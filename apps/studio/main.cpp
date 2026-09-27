// samaya Studio: the Windows app. A native window hosting Microsoft Edge WebView2; the page
// (apps/studio/ui) is the interface. Solves run out of process: each run starts samaya.exe (the
// CLI next to this program) with --json and --solution, streams its log to the page, then reads
// the solution file, re-checks it independently (studio_core) and sends the result. A solver
// crash therefore never takes the app down, and Stop always works.
//
//   samaya-studio.exe [--devtools] [--size 1600x1000]
//                     [--capture out.png [--run model.mps] [--tab plan] [--eval script] [--wait ms]]
//
// --capture renders the page to a PNG (CapturePreview) and exits; for automated checks.

#ifndef NOMINMAX
#define NOMINMAX
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <dwmapi.h>
#include <shellapi.h>
#include <shlobj.h>
#include <shlwapi.h>
#include <shobjidl.h>
#include <wrl.h>

#include <WebView2.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cwchar>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "studio_core.hpp"

using Microsoft::WRL::Callback;
using Microsoft::WRL::ComPtr;
namespace fs = std::filesystem;

namespace {

constexpr UINT kMsgPost = WM_APP + 1;  // lParam: std::string* JSON to post to the page
constexpr UINT_PTR kTimerCapture = 1;
constexpr UINT_PTR kTimerTab = 2;
constexpr wchar_t kHost[] = L"studio.samaya";
constexpr wchar_t kOrigin[] = L"https://studio.samaya/";  // the only page Studio shows

// UTF-8 <-> UTF-16 ------------------------------------------------------------------------------
std::wstring widen(const std::string& s) {
  if (s.empty()) return {};
  const int n = MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), nullptr, 0);
  std::wstring w(static_cast<std::size_t>(n), L'\0');
  MultiByteToWideChar(CP_UTF8, 0, s.data(), static_cast<int>(s.size()), w.data(), n);
  return w;
}
std::string narrow(const std::wstring& w) {
  if (w.empty()) return {};
  const int n = WideCharToMultiByte(CP_UTF8, 0, w.data(), static_cast<int>(w.size()), nullptr, 0, nullptr, nullptr);
  std::string s(static_cast<std::size_t>(n), '\0');
  WideCharToMultiByte(CP_UTF8, 0, w.data(), static_cast<int>(w.size()), s.data(), n, nullptr, nullptr);
  return s;
}
std::string read_file(const fs::path& p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

// State -----------------------------------------------------------------------------------------
struct Job {
  HANDLE process = nullptr;
  std::atomic<bool> cancelled{false};
  fs::path solution;
};

struct App {
  HWND hwnd = nullptr;
  ComPtr<ICoreWebView2Controller> controller;
  ComPtr<ICoreWebView2> webview;
  fs::path exe_dir;
  fs::path data_dir;  // %LOCALAPPDATA%\samaya\Studio
  HANDLE job_object = nullptr;
  bool dark = false;
  bool devtools = false;
  // Screenshot mode.
  std::wstring capture_png;
  std::string capture_run;
  std::string capture_tab;
  std::string capture_eval;  // script run on the page before the capture (for example studio.compare())
  int capture_wait = 700;    // milliseconds between that step and the capture
  int capture_stage = 0;     // 0: waiting for the first result; 1: past it
  int width = 1440, height = 920;  // window size in DIPs (--size)
  std::mutex mu;
  std::map<long long, std::shared_ptr<Job>> jobs;
  std::set<std::wstring> saved;  // files Studio wrote; the only ones the page may open or reveal
};
App g;

void post_to_page_now(const std::string& json) {
  if (g.webview) g.webview->PostWebMessageAsJson(widen(json).c_str());
}
// Thread-safe: hands the message to the UI thread.
void post_to_page(std::string json) {
  PostMessageW(g.hwnd, kMsgPost, 0, reinterpret_cast<LPARAM>(new std::string(std::move(json))));
}
void notice(const std::string& text) {
  post_to_page("{\"type\":\"notice\",\"message\":" + studio::json_string(text) + "}");
}

// The title bar in the page's colours (Windows 11; older versions ignore the attributes).
void style_title_bar() {
  const BOOL dark = g.dark ? TRUE : FALSE;
  DwmSetWindowAttribute(g.hwnd, 20 /* DWMWA_USE_IMMERSIVE_DARK_MODE */, &dark, sizeof dark);
  const COLORREF caption = g.dark ? RGB(0x23, 0x26, 0x2B) : RGB(0xFF, 0xFF, 0xFF);
  const COLORREF text = g.dark ? RGB(0xE3, 0xE6, 0xEA) : RGB(0x1D, 0x23, 0x2B);
  DwmSetWindowAttribute(g.hwnd, 35 /* DWMWA_CAPTION_COLOR */, &caption, sizeof caption);
  DwmSetWindowAttribute(g.hwnd, 36 /* DWMWA_TEXT_COLOR */, &text, sizeof text);
}

// Files -----------------------------------------------------------------------------------------
std::string files_json(const std::vector<std::wstring>& paths) {
  std::string out = "[";
  for (const auto& p : paths) {
    if (out.size() > 1) out += ",";
    out += "{\"path\":" + studio::json_string(narrow(p)) + ",\"name\":" +
           studio::json_string(narrow(fs::path(p).filename().wstring())) + "}";
  }
  return out + "]";
}

void pick_files() {
  ComPtr<IFileOpenDialog> dlg;
  if (FAILED(CoCreateInstance(CLSID_FileOpenDialog, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&dlg)))) return;
  const COMDLG_FILTERSPEC types[] = {{L"Optimization models (*.mps; *.qps)", L"*.mps;*.qps"}, {L"All files", L"*.*"}};
  dlg->SetFileTypes(2, types);
  dlg->SetTitle(L"Add models to solve");
  DWORD opts = 0;
  dlg->GetOptions(&opts);
  dlg->SetOptions(opts | FOS_ALLOWMULTISELECT | FOS_FILEMUSTEXIST | FOS_FORCEFILESYSTEM);
  if (FAILED(dlg->Show(g.hwnd))) return;
  ComPtr<IShellItemArray> items;
  if (FAILED(dlg->GetResults(&items))) return;
  DWORD count = 0;
  items->GetCount(&count);
  std::vector<std::wstring> paths;
  for (DWORD i = 0; i < count; ++i) {
    ComPtr<IShellItem> item;
    PWSTR path = nullptr;
    if (SUCCEEDED(items->GetItemAt(i, &item)) && SUCCEEDED(item->GetDisplayName(SIGDN_FILESYSPATH, &path))) {
      paths.emplace_back(path);
      CoTaskMemFree(path);
    }
  }
  post_to_page_now("{\"type\":\"picked\",\"files\":" + files_json(paths) + "}");
}

// Save dialog; returns the chosen path or empty.
std::wstring ask_save_path(const std::wstring& name, const std::wstring& ext) {
  ComPtr<IFileSaveDialog> dlg;
  if (FAILED(CoCreateInstance(CLSID_FileSaveDialog, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&dlg)))) return {};
  COMDLG_FILTERSPEC type = {L"File", L"*.*"};
  if (ext == L"xlsx") type = {L"Excel workbook (*.xlsx)", L"*.xlsx"};
  if (ext == L"zip") type = {L"Zip archive (*.zip)", L"*.zip"};
  if (ext == L"sol") type = {L"samaya solution (*.sol)", L"*.sol"};
  dlg->SetFileTypes(1, &type);
  dlg->SetDefaultExtension(ext.c_str());
  dlg->SetFileName(name.c_str());
  if (FAILED(dlg->Show(g.hwnd))) return {};
  ComPtr<IShellItem> item;
  PWSTR path = nullptr;
  if (FAILED(dlg->GetResult(&item)) || FAILED(item->GetDisplayName(SIGDN_FILESYSPATH, &path))) return {};
  std::wstring out(path);
  CoTaskMemFree(path);
  return out;
}

void samples_json(std::string& out) {
  struct Known {
    const char* file;
    const char* label;
  };
  const Known known[] = {{"mrpl_crude_small.mps", "Crude receipt schedule"},
                         {"mrpl_plan_small.mps", "Refinery plan"},
                         {"mrpl_utility_small.mps", "Steam and power"},
                         {"mrpl_crude_medium.mps", "Crude schedule (21 days)"},
                         {"mrpl_plan_medium.mps", "Refinery plan (13 weeks)"}};
  out = "[";
  for (const Known& k : known) {
    const fs::path p = g.exe_dir / "samples" / k.file;
    if (!fs::exists(p)) continue;
    if (out.size() > 1) out += ",";
    out += "{\"path\":" + studio::json_string(narrow(p.wstring())) + ",\"name\":" + studio::json_string(k.file) +
           ",\"label\":" + studio::json_string(k.label) + "}";
  }
  out += "]";
}

// Solving ---------------------------------------------------------------------------------------
std::wstring quote(const std::wstring& s) { return L"\"" + s + L"\""; }

void solve(long long id, const std::string& model, const studio::Json& options) {
  auto job = std::make_shared<Job>();
  const fs::path run_dir = g.data_dir / "runs" / std::to_wstring(GetCurrentProcessId()) / std::to_wstring(id);
  std::error_code ec;
  fs::create_directories(run_dir, ec);
  job->solution = run_dir / L"solution.sol";
  const fs::path exe = g.exe_dir / L"samaya.exe";
  if (!fs::exists(exe)) {
    post_to_page("{\"type\":\"error\",\"id\":" + std::to_string(id) + ",\"message\":\"samaya.exe was not found next to Studio.\"}");
    return;
  }
  const std::string method = options.str("lpMethod", "auto");
  std::wstring cmd = quote(exe.wstring()) + L" --json --solution " + quote(job->solution.wstring());
  cmd += L" --time-limit " + std::to_wstring(options.num("timeLimit", 60));
  const int threads = static_cast<int>(options.num("threads", 0));
  if (threads > 0) cmd += L" --threads " + std::to_wstring(threads);
  cmd += L" --mip-gap " + widen(studio::json_number(options.num("mipGap", 1e-4)));
  cmd += L" --lp-method " + widen(method);
  if (!options.flag("presolve", true)) cmd += L" --no-presolve";
  cmd += L" " + quote(widen(model));

  SECURITY_ATTRIBUTES sa{sizeof sa, nullptr, TRUE};
  HANDLE read_end = nullptr, write_end = nullptr;
  if (!CreatePipe(&read_end, &write_end, &sa, 0)) {  // without an answer the page would wait forever
    post_to_page("{\"type\":\"error\",\"id\":" + std::to_string(id) + ",\"message\":\"Could not start samaya.exe.\"}");
    return;
  }
  SetHandleInformation(read_end, HANDLE_FLAG_INHERIT, 0);
  STARTUPINFOW si{};
  si.cb = sizeof si;
  si.dwFlags = STARTF_USESTDHANDLES;
  si.hStdOutput = write_end;
  si.hStdError = write_end;
  si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
  PROCESS_INFORMATION pi{};
  std::vector<wchar_t> buf(cmd.begin(), cmd.end());
  buf.push_back(L'\0');
  const BOOL started = CreateProcessW(nullptr, buf.data(), nullptr, nullptr, TRUE, CREATE_NO_WINDOW | CREATE_SUSPENDED,
                                      nullptr, g.exe_dir.c_str(), &si, &pi);
  CloseHandle(write_end);
  if (!started) {
    CloseHandle(read_end);
    post_to_page("{\"type\":\"error\",\"id\":" + std::to_string(id) + ",\"message\":\"Could not start samaya.exe.\"}");
    return;
  }
  AssignProcessToJobObject(g.job_object, pi.hProcess);  // dies with Studio
  ResumeThread(pi.hThread);
  CloseHandle(pi.hThread);
  job->process = pi.hProcess;
  {
    std::lock_guard<std::mutex> lock(g.mu);
    g.jobs[id] = job;
  }

  std::thread([id, job, read_end, model]() {
    const auto start = std::chrono::steady_clock::now();
    std::string pending, all;
    char chunk[4096];
    DWORD got = 0;
    while (ReadFile(read_end, chunk, sizeof chunk, &got, nullptr) && got > 0) {
      pending.append(chunk, got);
      std::size_t nl;
      while ((nl = pending.find('\n')) != std::string::npos) {
        std::string line = pending.substr(0, nl);
        pending.erase(0, nl + 1);
        if (!line.empty() && line.back() == '\r') line.pop_back();
        all += line + "\n";
        if (!line.empty() && line[0] != '{') {  // the JSON summary is not log text
          post_to_page("{\"type\":\"log\",\"id\":" + std::to_string(id) + ",\"line\":" + studio::json_string(line) + "}");
        }
      }
    }
    if (!pending.empty()) all += pending;
    CloseHandle(read_end);
    WaitForSingleObject(job->process, INFINITE);
    DWORD code = 0;
    GetExitCodeProcess(job->process, &code);
    {
      // Under the lock, so that a late Stop never terminates a closed (possibly reused) handle.
      std::lock_guard<std::mutex> lock(g.mu);
      CloseHandle(job->process);
      job->process = nullptr;
    }
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    const std::string sid = std::to_string(id);
    if (job->cancelled) {
      post_to_page("{\"type\":\"cancelled\",\"id\":" + sid + "}");
      return;
    }
    // The last JSON line of the output is the run summary; the rest is the log.
    std::string summary, log;
    std::istringstream lines(all);
    for (std::string line; std::getline(lines, line);) {
      if (!line.empty() && line[0] == '{') summary = line;
      else log += line + "\n";
    }
    if (summary.empty()) {
      std::string msg = log.empty() ? "samaya stopped with exit code " + std::to_string(code) : log;
      post_to_page("{\"type\":\"error\",\"id\":" + sid + ",\"message\":" + studio::json_string(msg) + "}");
      return;
    }
    std::string result = "{\"type\":\"result\",\"id\":" + sid + ",\"summary\":" + summary;
    result += ",\"elapsed\":" + studio::json_number(elapsed);
    result += ",\"path\":" + studio::json_string(model) + ",\"log\":" + studio::json_string(log);
    // The solution and the independent re-check.
    studio::SolutionFile sol;
    try {
      sol = studio::read_solution(read_file(job->solution));
    } catch (const std::exception&) {
    }
    studio::Recheck chk;
    try {
      chk = studio::recheck(read_file(fs::path(widen(model))), sol);
    } catch (const std::exception& e) {
      chk.reason = e.what();
    }
    result += ",\"sense\":" + studio::json_string(chk.maximize ? "maximize" : "minimize");
    result += ",\"recheck\":{\"available\":" + std::string(chk.available ? "true" : "false") +
              ",\"reason\":" + studio::json_string(chk.reason) + ",\"rows\":" + studio::json_number(chk.row_violation) +
              ",\"rowName\":" + studio::json_string(chk.row_name) + ",\"bounds\":" + studio::json_number(chk.bound_violation) +
              ",\"integrality\":" + studio::json_number(chk.integrality) + ",\"objective\":" + studio::json_number(chk.objective) +
              ",\"objectiveDiff\":" + studio::json_number(chk.quadratic ? 0.0 : chk.objective_diff) +
              ",\"quadratic\":" + std::string(chk.quadratic ? "true" : "false") +
              ",\"missing\":" + std::to_string(chk.missing) + "}";
    // At most this many entries go to the page; the solution file keeps every value.
    constexpr std::size_t kMaxEntries = 200000;
    result += ",\"cols\":[";
    for (std::size_t j = 0; j < sol.columns.size() && j < kMaxEntries; ++j) {
      if (j) result += ",";
      result += "[" + studio::json_string(sol.columns[j].first) + "," + studio::json_number(sol.columns[j].second) + "]";
    }
    result += "],\"rowsData\":[";
    for (std::size_t i = 0; i < sol.rows.size() && i < kMaxEntries; ++i) {
      if (i) result += ",";
      const auto& r = sol.rows[i];
      result += "[" + studio::json_string(r.name) + "," + studio::json_number(r.activity) + "," + studio::json_number(r.dual) + "]";
    }
    result += "],\"truncated\":" + std::string(sol.columns.size() > kMaxEntries || sol.rows.size() > kMaxEntries ? "true" : "false");
    result += ",\"solutionPath\":" + studio::json_string(narrow(job->solution.wstring())) + "}";
    post_to_page(std::move(result));
  }).detach();
}

void cancel(long long id) {
  std::lock_guard<std::mutex> lock(g.mu);
  auto it = g.jobs.find(id);
  if (it == g.jobs.end() || !it->second->process) return;  // unknown, or already finished
  it->second->cancelled = true;
  TerminateProcess(it->second->process, 1);
}

// Messages from the page ---------------------------------------------------------------------
void on_message(ICoreWebView2WebMessageReceivedEventArgs* args) {
  LPWSTR source = nullptr;
  if (FAILED(args->get_Source(&source))) return;
  const bool ours = std::wstring_view(source).starts_with(kOrigin);
  CoTaskMemFree(source);
  if (!ours) return;
  LPWSTR raw = nullptr;
  if (FAILED(args->get_WebMessageAsJson(&raw))) return;
  const std::string text = narrow(raw);
  CoTaskMemFree(raw);
  studio::Json m;
  try {
    m = studio::parse_json(text);
  } catch (const std::exception&) {
    return;
  }
  const std::string type = m.str("type");
  const long long id = static_cast<long long>(m.num("id", 0));
  if (type == "ready") {
    std::string samples;
    samples_json(samples);
    post_to_page_now("{\"type\":\"init\",\"version\":" + studio::json_string(SAMAYA_VERSION_STRING) +
                     ",\"cores\":" + std::to_string(std::thread::hardware_concurrency()) +
                     ",\"gpu\":false,\"samples\":" + samples + "}");
    if (!g.capture_run.empty()) {
      post_to_page_now("{\"type\":\"picked\",\"files\":" + files_json({widen(g.capture_run)}) + "}");
    } else if (!g.capture_png.empty()) {
      g.capture_stage = 1;
      SetTimer(g.hwnd, g.capture_tab.empty() && g.capture_eval.empty() ? kTimerCapture : kTimerTab, 1500, nullptr);
    }
  } else if (type == "theme") {
    // The page chooses the theme (light by default); the title bar and the background follow it.
    g.dark = m.flag("dark", false);
    style_title_bar();
    ComPtr<ICoreWebView2Controller2> c2;
    if (g.controller && SUCCEEDED(g.controller.As(&c2))) {
      c2->put_DefaultBackgroundColor(g.dark ? COREWEBVIEW2_COLOR{255, 0x1B, 0x1D, 0x21} : COREWEBVIEW2_COLOR{255, 0xF3, 0xF4, 0xF6});
    }
  } else if (type == "pick") {
    pick_files();
  } else if (type == "dropped") {
    // The page passes the dropped File objects; WebView2 gives their real paths.
    ComPtr<ICoreWebView2WebMessageReceivedEventArgs2> args2;
    ComPtr<ICoreWebView2ObjectCollectionView> objects;
    UINT32 count = 0;
    std::vector<std::wstring> paths;
    if (SUCCEEDED(args->QueryInterface(IID_PPV_ARGS(&args2))) && SUCCEEDED(args2->get_AdditionalObjects(&objects)) &&
        SUCCEEDED(objects->get_Count(&count))) {
      for (UINT32 i = 0; i < count; ++i) {
        ComPtr<IUnknown> obj;
        ComPtr<ICoreWebView2File> file;
        LPWSTR path = nullptr;
        if (SUCCEEDED(objects->GetValueAtIndex(i, &obj)) && obj && SUCCEEDED(obj.As(&file)) &&
            SUCCEEDED(file->get_Path(&path))) {
          paths.emplace_back(path);
          CoTaskMemFree(path);
        }
      }
    }
    post_to_page_now("{\"type\":\"picked\",\"files\":" + files_json(paths) + "}");
  } else if (type == "solve") {
    const studio::Json* opts = m.get("options");
    solve(id, m.str("path"), opts ? *opts : studio::Json());
  } else if (type == "cancel") {
    cancel(id);
  } else if (type == "save") {
    const std::wstring name = widen(m.str("name", "samaya.xlsx"));
    const std::wstring ext = m.str("kind") == "zip" ? L"zip" : L"xlsx";
    const std::wstring path = ask_save_path(name, ext);
    if (path.empty()) return;
    const std::string bytes = studio::base64_decode(m.str("base64"));
    std::ofstream out(fs::path(path), std::ios::binary);
    out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!out) {
      notice("Could not write the file.");
      return;
    }
    g.saved.insert(path);
    post_to_page_now("{\"type\":\"saved\",\"path\":" + studio::json_string(narrow(path)) + "}");
  } else if (type == "saveCopy") {
    fs::path from;
    {
      std::lock_guard<std::mutex> lock(g.mu);
      auto it = g.jobs.find(id);
      if (it != g.jobs.end()) from = it->second->solution;
    }
    if (from.empty() || !fs::exists(from)) {
      notice("This run has no solution file.");
      return;
    }
    const std::wstring path = ask_save_path(widen(m.str("name", "solution.sol")), L"sol");
    if (path.empty()) return;
    std::error_code ec;
    fs::copy_file(from, fs::path(path), fs::copy_options::overwrite_existing, ec);
    if (ec) {
      notice("Could not save the solution file.");
      return;
    }
    g.saved.insert(path);
    post_to_page_now("{\"type\":\"saved\",\"path\":" + studio::json_string(narrow(path)) + "}");
  } else if (type == "open" || type == "reveal") {
    // Only a file Studio itself saved: the page never gets to open an arbitrary path.
    const std::wstring path = widen(m.str("path"));
    if (!g.saved.contains(path)) return;
    if (type == "open") {
      ShellExecuteW(g.hwnd, L"open", path.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
      return;
    }
    const std::wstring arg = L"/select,\"" + path + L"\"";
    ShellExecuteW(g.hwnd, nullptr, L"explorer.exe", arg.c_str(), nullptr, SW_SHOWNORMAL);
  }
}

// Screenshot mode: CapturePreview writes the page as PNG, then the app exits.
void capture() {
  ComPtr<IStream> stream;
  if (FAILED(SHCreateStreamOnFileEx(g.capture_png.c_str(), STGM_CREATE | STGM_WRITE, FILE_ATTRIBUTE_NORMAL, TRUE,
                                    nullptr, &stream))) {
    PostQuitMessage(2);
    return;
  }
  IStream* raw = stream.Detach();
  g.webview->CapturePreview(COREWEBVIEW2_CAPTURE_PREVIEW_IMAGE_FORMAT_PNG, raw,
                            Callback<ICoreWebView2CapturePreviewCompletedHandler>([raw](HRESULT hr) -> HRESULT {
                              raw->Release();
                              PostQuitMessage(SUCCEEDED(hr) ? 0 : 3);
                              return S_OK;
                            }).Get());
}

// WebView2 setup -------------------------------------------------------------------------------
void create_webview() {
  const fs::path user_data = g.data_dir / L"WebView2";
  const HRESULT hr = CreateCoreWebView2EnvironmentWithOptions(
      nullptr, user_data.c_str(), nullptr,
      Callback<ICoreWebView2CreateCoreWebView2EnvironmentCompletedHandler>(
          [](HRESULT result, ICoreWebView2Environment* env) -> HRESULT {
            if (FAILED(result)) return result;
            env->CreateCoreWebView2Controller(
                g.hwnd, Callback<ICoreWebView2CreateCoreWebView2ControllerCompletedHandler>(
                            [](HRESULT result2, ICoreWebView2Controller* controller) -> HRESULT {
                              if (FAILED(result2) || !controller) return result2;
                              g.controller = controller;
                              g.controller->get_CoreWebView2(&g.webview);
                              ComPtr<ICoreWebView2Controller2> c2;
                              if (SUCCEEDED(g.controller.As(&c2))) {
                                COREWEBVIEW2_COLOR bg = g.dark ? COREWEBVIEW2_COLOR{255, 0x1B, 0x1D, 0x21}
                                                               : COREWEBVIEW2_COLOR{255, 0xF3, 0xF4, 0xF6};
                                c2->put_DefaultBackgroundColor(bg);
                              }
                              ComPtr<ICoreWebView2Settings> settings;
                              g.webview->get_Settings(&settings);
                              settings->put_IsStatusBarEnabled(FALSE);
                              settings->put_IsZoomControlEnabled(FALSE);
                              settings->put_AreDevToolsEnabled(g.devtools ? TRUE : FALSE);
                              settings->put_AreDefaultContextMenusEnabled(g.devtools ? TRUE : FALSE);
                              ComPtr<ICoreWebView2Settings3> s3;
                              if (SUCCEEDED(settings.As(&s3))) s3->put_AreBrowserAcceleratorKeysEnabled(g.devtools ? TRUE : FALSE);
                              RECT rc;
                              GetClientRect(g.hwnd, &rc);
                              g.controller->put_Bounds(rc);
                              ComPtr<ICoreWebView2_3> w3;
                              if (SUCCEEDED(g.webview.As(&w3))) {
                                w3->SetVirtualHostNameToFolderMapping(kHost, (g.exe_dir / L"ui").c_str(),
                                                                      COREWEBVIEW2_HOST_RESOURCE_ACCESS_KIND_DENY_CORS);
                              }
                              EventRegistrationToken token;
                              // Studio shows only its own page: no navigation elsewhere (a dropped
                              // file or a link), and no new windows.
                              g.webview->add_NavigationStarting(
                                  Callback<ICoreWebView2NavigationStartingEventHandler>(
                                      [](ICoreWebView2*, ICoreWebView2NavigationStartingEventArgs* a) -> HRESULT {
                                        LPWSTR uri = nullptr;
                                        if (SUCCEEDED(a->get_Uri(&uri))) {
                                          if (!std::wstring_view(uri).starts_with(kOrigin)) a->put_Cancel(TRUE);
                                          CoTaskMemFree(uri);
                                        }
                                        return S_OK;
                                      }).Get(),
                                  &token);
                              g.webview->add_NewWindowRequested(
                                  Callback<ICoreWebView2NewWindowRequestedEventHandler>(
                                      [](ICoreWebView2*, ICoreWebView2NewWindowRequestedEventArgs* a) -> HRESULT {
                                        a->put_Handled(TRUE);
                                        return S_OK;
                                      }).Get(),
                                  &token);
                              g.webview->add_WebMessageReceived(
                                  Callback<ICoreWebView2WebMessageReceivedEventHandler>(
                                      [](ICoreWebView2*, ICoreWebView2WebMessageReceivedEventArgs* args) -> HRESULT {
                                        on_message(args);
                                        return S_OK;
                                      }).Get(),
                                  &token);
                              g.webview->Navigate((std::wstring(L"https://") + kHost + L"/index.html").c_str());
                              return S_OK;
                            }).Get());
            return S_OK;
          }).Get());
  if (FAILED(hr)) {
    const int answer = MessageBoxW(g.hwnd,
                                   L"samaya Studio needs the Microsoft Edge WebView2 Runtime, which is part of Windows 11 "
                                   L"and current Windows 10.\n\nRun setup.bat, or choose OK to open Microsoft's download page.",
                                   L"samaya Studio", MB_OKCANCEL | MB_ICONINFORMATION);
    if (answer == IDOK) {
      ShellExecuteW(nullptr, L"open", L"https://go.microsoft.com/fwlink/p/?LinkId=2124703", nullptr, nullptr, SW_SHOWNORMAL);
    }
    PostQuitMessage(1);
  }
}

LRESULT CALLBACK window_proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
  switch (msg) {
    case WM_SIZE:
      if (g.controller) {
        RECT rc;
        GetClientRect(hwnd, &rc);
        g.controller->put_Bounds(rc);
      }
      return 0;
    case kMsgPost: {
      std::unique_ptr<std::string> json(reinterpret_cast<std::string*>(lp));
      post_to_page_now(*json);
      // Screenshot mode: once the first result is on the page, open the tab and run the script,
      // then capture.
      if (!g.capture_png.empty() && g.capture_stage == 0 && json->rfind("{\"type\":\"result\"", 0) == 0) {
        g.capture_stage = 1;
        SetTimer(hwnd, g.capture_tab.empty() && g.capture_eval.empty() ? kTimerCapture : kTimerTab, 900, nullptr);
      }
      return 0;
    }
    case WM_TIMER:
      KillTimer(hwnd, wp);
      if (wp == kTimerTab && g.webview) {
        if (!g.capture_tab.empty()) {
          g.webview->ExecuteScript(widen("studio.select(" + studio::json_string(g.capture_tab) + ")").c_str(), nullptr);
        }
        if (!g.capture_eval.empty()) g.webview->ExecuteScript(widen(g.capture_eval).c_str(), nullptr);
        SetTimer(hwnd, kTimerCapture, static_cast<UINT>(g.capture_wait), nullptr);
      } else if (wp == kTimerCapture && g.webview) {
        capture();
      }
      return 0;
    case WM_DESTROY:
      PostQuitMessage(0);
      return 0;
    default:
      return DefWindowProcW(hwnd, msg, wp, lp);
  }
}

}  // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show) {
  SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
  if (FAILED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED))) return 1;

  wchar_t exe[MAX_PATH];
  GetModuleFileNameW(nullptr, exe, MAX_PATH);
  g.exe_dir = fs::path(exe).parent_path();
  PWSTR local = nullptr;
  if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_LocalAppData, 0, nullptr, &local))) {
    g.data_dir = fs::path(local) / L"samaya" / L"Studio";
    CoTaskMemFree(local);
  } else {
    g.data_dir = fs::temp_directory_path() / L"samaya-studio";
  }
  std::error_code ec;
  fs::create_directories(g.data_dir, ec);

  int argc = 0;
  LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
  for (int i = 1; i < argc; ++i) {
    const std::wstring a = argv[i];
    if (a == L"--devtools") g.devtools = true;
    if (a == L"--capture" && i + 1 < argc) g.capture_png = argv[++i];
    if (a == L"--run" && i + 1 < argc) g.capture_run = narrow(argv[++i]);
    if (a == L"--tab" && i + 1 < argc) g.capture_tab = narrow(argv[++i]);
    if (a == L"--eval" && i + 1 < argc) g.capture_eval = narrow(argv[++i]);
    if (a == L"--wait" && i + 1 < argc) g.capture_wait = std::clamp(_wtoi(argv[++i]), 0, 60000);
    if (a == L"--size" && i + 1 < argc) {
      int sw = 0, sh = 0;
      if (swscanf_s(argv[++i], L"%dx%d", &sw, &sh) == 2 && sw >= 640 && sh >= 480) {
        g.width = sw;
        g.height = sh;
      }
    }
  }
  LocalFree(argv);

  // Child solver processes die with Studio.
  g.job_object = CreateJobObjectW(nullptr, nullptr);
  JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
  limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
  SetInformationJobObject(g.job_object, JobObjectExtendedLimitInformation, &limits, sizeof limits);

  g.dark = false;  // light until the page asks for dark (its choice is kept between starts)
  WNDCLASSEXW wc{};
  wc.cbSize = sizeof wc;
  wc.lpfnWndProc = window_proc;
  wc.hInstance = instance;
  wc.hIcon = LoadIconW(instance, L"IDI_APP");
  wc.hIconSm = wc.hIcon;
  wc.hCursor = LoadCursorW(nullptr, IDC_ARROW);
  wc.hbrBackground = CreateSolidBrush(RGB(0xF3, 0xF4, 0xF6));
  wc.lpszClassName = L"SamayaStudio";
  RegisterClassExW(&wc);

  const UINT dpi = GetDpiForSystem();
  const int w = MulDiv(g.width, static_cast<int>(dpi), 96), h = MulDiv(g.height, static_cast<int>(dpi), 96);
  RECT work;
  SystemParametersInfoW(SPI_GETWORKAREA, 0, &work, 0);
  const int x = work.left + std::max(0, static_cast<int>((work.right - work.left - w) / 2));
  const int y = work.top + std::max(0, static_cast<int>((work.bottom - work.top - h) / 2));
  g.hwnd = CreateWindowExW(0, wc.lpszClassName, L"samaya Studio", WS_OVERLAPPEDWINDOW, x, y,
                           std::min(w, static_cast<int>(work.right - work.left)),
                           std::min(h, static_cast<int>(work.bottom - work.top)), nullptr, nullptr, instance, nullptr);
  if (!g.hwnd) return 1;
  style_title_bar();
  ShowWindow(g.hwnd, show);
  UpdateWindow(g.hwnd);
  create_webview();

  MSG msg;
  while (GetMessageW(&msg, nullptr, 0, 0) > 0) {
    TranslateMessage(&msg);
    DispatchMessageW(&msg);
  }
  g.controller.Reset();
  g.webview.Reset();
  if (g.job_object) CloseHandle(g.job_object);
  CoUninitialize();
  return static_cast<int>(msg.wParam);
}
