# Plan: samaya Studio for refinery planners, calendar dates, and the demo video (28 Sep 2026, night)

The owner asked, before the SIH demo on the morning of 28 September:
- a GUI shaped by how refinery planners and schedulers actually work, compared with the leading
  products;
- no flashy "AI" styling: professional and premium;
- calendar dates instead of "day 1, day 2";
- then a full re-check of the project for bugs;
- and finally Studio in the demo video.

## Status (28 Sep 2026, 02:30 IST)

**Done:**
- **P1–P4:**
  - the new design: a toolbar, a list of cases, a status bar, light by default and dark on
    request;
  - the start date and calendar axes;
  - real Excel dates, with the header row frozen;
  - Compare cases;
  - the berth Gantt chart, tank stocks, CDU charge, weekly plan and hourly unit commitment.
- **P5:**
  - every page checked in the real app (screenshots of the three samples, Compare, dark theme,
    a malformed file);
  - Excel opened the exported workbook and shows the dates, for example `Thu 01-Oct-2026`;
  - all six presets passed, and the package was rebuilt and installed.
- **P6 is outside the repo:** the video's demo scene now shows the installed Studio (see the
  report to the owner).

**Not done:**
- Tank capacity lines on the stock chart: the page does not get column bounds from the solution
  file.
- Utilization against capacity in the plan: the same reason.



**How MRPL uses this kind of software.**
- MRPL is a 15 MMTPA refinery that runs crudes from 18 to 46 API
  ([MRPL refining](https://mrpl.co.in/en/Content/Refining)).
- Its planners run LP models for the crude slate, unit loadings and product economics. Its
  schedulers turn the plan into crude receipts, tank transfers and CDU feed.
- The problem statement (SIH26119) asks for an indigenous LP, MILP and QP solver in place of
  Xpress and CPLEX.
- So Studio's users are **planners and schedulers who today use Aspen PIMS and a scheduler**,
  with Excel as the common language.

**What the leading products show.**

| Product | Users | What its screens are built around |
|---|---|---|
| Aspen PIMS / PIMS Platinum ([brochure](https://www.aspentech.com/en/resources/brochure/aspen-pims-family)) | Planners, at more than 400 refineries | **Cases**: a base case and alternatives run side by side ("evaluate multiple cases at one time"); crude slate; unit utilization; marginal values; Excel reports |
| Aspen Petroleum Scheduler ([product](https://www.aspentech.com/en/products/msc/aspen-petroleum-scheduler)) | Schedulers, at more than 250 sites | **Gantt charts** of crude receipt, tank transfer and CDU feed events, on a calendar; flowsheets |
| Haverly H/Sched ([brochure](https://www.haverly.com/images/downloads/products/about%20hsched%202021.pdf)) | Schedulers | "The window a user works from most is the one containing the schedule **Gantt charts**"; **tank graphs** by quantity and quality; results sent to **Excel** reports and work orders |
| FICO Xpress Insight ([overview](https://www.fico.com/en/products/fico-xpress-insight)), GAMS MIRO ([docs](https://www.gams.com/miro/start.html)) | Business users of optimization models | **Scenarios**: inputs and outputs per run, and comparison in split, tab and pivot views |
| ISA-101 high-performance HMI ([summary](https://processcontrolguide.com/isa-101-hmi-design/)) | Refinery control rooms | Mostly neutral grey; **colour only for meaning** (normal, advisory, alarm); no decorative graphics; a trend within one click |

**What this means for Studio:**
1. **Calendar time.** Every period is a date:
   - the crude schedule is days;
   - the plan is weeks (shown as "week of …");
   - the utility case is hours on dates.
   - A **Start date** setting anchors them. The Gantt axis shows dates, weekdays and month
     boundaries, with weekends shaded.
2. **Cases, not "runs".** A run is a case, with its settings shown next to its result. **Compare
   cases** puts two or more runs of the same model side by side (objective, status, time,
   settings, and the largest differences in the plan), as in PIMS and MIRO. Back / Run again
   (already built) creates the alternatives.
3. **Schedule first for scheduling models.** For the crude case, the Plan tab opens with the berth
   Gantt, then tank stocks against capacity, then CDU charge. Every chart has its table, and
   Excel gets real date cells.
4. **Sober, dense and neutral**, following ISA-101 and Windows (Fluent) conventions:
   - Segoe UI at 13 px, with tabular figures;
   - a grey and white base, and one dark-blue accent for actions;
   - green, amber and red only for status;
   - no gradients, glows, pill-shaped badges, emoji check marks or large hero text;
   - a toolbar with the main commands (Open, Run, Stop, Export), and a status bar at the bottom;
   - light theme by default, as office software.
5. **Verification stays visible but factual:** a "Verified" status with the numbers behind it, in
   the same sober style.

## 2. Work, in order

- **P1, design system (`apps/studio/ui/app.css`, `index.html`):**
  - neutral tokens, the type scale, and the toolbar, sidebar, tabs, table, card and status bar
    styles;
  - remove the gradients, the hero, the rounded pills and the emoji;
  - light theme by default, with dark mode only when Windows is dark.
- **P2, calendar (`app.js`, `charts.js`):**
  - a Start date setting (default: tomorrow), kept per case;
  - a period model per case (day, week or hour);
  - formatting in Indian style (`Mon 28 Sep 2026`, `28-Sep-2026`);
  - date axes with weekend shading and month labels;
  - tables and Excel with dates (real Excel date cells in `xlsx.js`).
- **P3, cases:**
  - rename runs to cases in the UI;
  - a settings summary on each case;
  - **Compare** (select two or more cases of the same model): KPIs, settings, and the largest
    differences in variables, in the app and in Excel.
- **P4, the scheduler views:**
  - a berth Gantt with a bar per cargo;
  - tank stocks against their capacity lines;
  - CDU charge stacked by crude;
  - for the plan case: crude slate by week, unit utilization against capacity, and marginal
    values;
  - for the utility case: unit on/off by hour.
- **P5, verification (no shortcuts):**
  - screenshots of every page in both themes through Studio's `--capture` mode, looked at one by
    one;
  - the three samples, a malformed file, an infeasible model and a QPS model in the real app;
  - the Excel workbook opened in Excel;
  - Linux debug, release and asan, and Windows debug, release and asan, all green with zero
    warnings;
  - the package rebuilt and installed with `setup.bat`;
  - one commit with the exact test output. Nothing is pushed.
- **P6, the video (`~/video-work`):**
  - replace the demo scene's browser report with the real samaya Studio: a Windows 11 frame
    around Studio's own screenshots;
  - dropping the crude case, the result, the calendar berth schedule, the verification, and Excel;
  - new narration for that scene only (the same Indian English voice and pace), new captions, a
    new music bed for the new length;
  - render, check stills of every scene, write the `.srt`, and replace the Desktop video (the
    current one is kept as "…(report demo, old)").

## 3. Not in tonight's scope

Editing model data in the app (Studio solves MPS files; building models stays in `cases/mrpl.py`
and the planners' own tools), multi-user work, and a database of cases.
