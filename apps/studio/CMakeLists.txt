# samaya Studio: the Windows desktop app. A Win32 window hosting Microsoft Edge WebView2 (part of
# Windows 11 and current Windows 10); the interface is the page in ui/. Solves run in samaya.exe.

# The WebView2 SDK (headers and the static loader), pinned by version and SHA-256. For offline
# builds pass -DSAMAYA_WEBVIEW2_NUPKG=<path to the .nupkg>.
set(SAMAYA_WEBVIEW2_VERSION 1.0.4191.47)
set(SAMAYA_WEBVIEW2_SHA256 f492bbf547d0da329553b6727435b677579b1e9f91cc9e4a1ad029366d5f23d0)
set(SAMAYA_WEBVIEW2_NUPKG "" CACHE FILEPATH "Local copy of the Microsoft.Web.WebView2 .nupkg")
set(wv2_dir ${CMAKE_BINARY_DIR}/_deps/webview2-${SAMAYA_WEBVIEW2_VERSION})
if(NOT EXISTS ${wv2_dir}/build/native/include/WebView2.h)
  set(wv2_pkg ${CMAKE_BINARY_DIR}/_deps/webview2-${SAMAYA_WEBVIEW2_VERSION}.nupkg)
  if(SAMAYA_WEBVIEW2_NUPKG)
    configure_file(${SAMAYA_WEBVIEW2_NUPKG} ${wv2_pkg} COPYONLY)
  else()
    file(DOWNLOAD
      https://api.nuget.org/v3-flatcontainer/microsoft.web.webview2/${SAMAYA_WEBVIEW2_VERSION}/microsoft.web.webview2.${SAMAYA_WEBVIEW2_VERSION}.nupkg
      ${wv2_pkg} STATUS wv2_status)
    list(GET wv2_status 0 wv2_code)
    if(NOT wv2_code EQUAL 0)
      message(FATAL_ERROR "Could not download the WebView2 SDK (${wv2_status}); pass SAMAYA_WEBVIEW2_NUPKG.")
    endif()
  endif()
  file(SHA256 ${wv2_pkg} wv2_hash)
  if(NOT wv2_hash STREQUAL SAMAYA_WEBVIEW2_SHA256)
    message(FATAL_ERROR "WebView2 SDK checksum mismatch: ${wv2_hash}")
  endif()
  file(ARCHIVE_EXTRACT INPUT ${wv2_pkg} DESTINATION ${wv2_dir} PATTERNS "build/native/*")
endif()

add_executable(samaya_studio WIN32 main.cpp studio_core.cpp studio.rc studio.manifest)
set_target_properties(samaya_studio PROPERTIES OUTPUT_NAME samaya-studio)
target_include_directories(samaya_studio PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${wv2_dir}/build/native/include)
target_compile_definitions(samaya_studio PRIVATE UNICODE _UNICODE SAMAYA_VERSION_STRING="${PROJECT_VERSION}")
target_link_libraries(samaya_studio PRIVATE samaya_options
  ${wv2_dir}/build/native/x64/WebView2LoaderStatic.lib dwmapi shlwapi version)
add_dependencies(samaya_studio samaya_cli)
# Next to the app: the solver it runs and the page.
add_custom_command(TARGET samaya_studio POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:samaya_cli> $<TARGET_FILE_DIR:samaya_studio>/samaya.exe
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/ui $<TARGET_FILE_DIR:samaya_studio>/ui)
install(TARGETS samaya_studio RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY ui DESTINATION ${CMAKE_INSTALL_BINDIR})
