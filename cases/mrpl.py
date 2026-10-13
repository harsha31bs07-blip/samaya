#!/usr/bin/env python3
"""MRPL-style refinery case studies: model generators and a readable plan report.

The figures are public, order-of-magnitude values for a 15 MMTPA coastal refinery (about
45 kt of crude a day), not MRPL data. Quantities are in kilotonnes (kt), money in lakh rupees
(1 lakh = 100,000), time in days (planning, scheduling) or hours (utilities).

Families:
  plan      multi-period refinery planning (LP): crude slate from a set of assayed crudes;
            CDU with crude-dependent straight-run yields and a blend sulfur limit; VDU;
            FCC, hydrocracker, naphtha reformer and diesel hydrotreater (DHDS); blending to BS-VI
            specifications (diesel sulfur <= 10 ppm, gasoline RON >= 91 and sulfur <= 10 ppm);
            product tanks; domestic demand and export at a discount.
  crude     crude receipt scheduling (MILP): cargoes with arrival windows at one single-point
            mooring (one cargo a day), optional spot cargoes, crude tanks, a CDU run between
            minimum and maximum throughput with a sulfur limit on the blend, demurrage.
  utility   captive power and steam unit commitment (MILP): boilers, a gas turbine with a heat
            recovery steam generator, a back-pressure steam turbine, a letdown valve and grid
            import with time-of-day tariffs; on/off and startup decisions, minimum up time,
            minimum and maximum loads, ramp limits, over 24 or 48 hours.

Usage:
  cases/mrpl.py generate --out cases/instances            # all families, three sizes each
  cases/mrpl.py report   MODEL.mps SOLUTION.txt            # readable plan from a samaya solution
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "bench"))
from generate_lps import MpsWriter  # noqa: E402

# name, API gravity, sulfur (wt %), cost (lakh Rs per kt, about $560-640 per tonne)
CRUDES = [
    ("arab_light", 33.0, 1.8, 4850.0),
    ("arab_heavy", 28.0, 2.9, 4600.0),
    ("basrah_medium", 29.0, 2.9, 4580.0),
    ("murban", 40.0, 0.8, 5150.0),
    ("upper_zakum", 34.0, 1.9, 4880.0),
    ("kuwait_export", 31.0, 2.6, 4700.0),
    ("bonny_light", 35.0, 0.14, 5250.0),
    ("das_blend", 39.0, 1.1, 5050.0),
]

# Straight-run cuts. Yields interpolate between a heavy (API 28) and a light (API 40) crude.
CUTS = ["lpg", "naphtha", "kero", "gasoil", "vgo", "residue"]
YIELD_HEAVY = [0.010, 0.120, 0.110, 0.220, 0.250, 0.290]
YIELD_LIGHT = [0.025, 0.230, 0.150, 0.250, 0.200, 0.145]

# Products: name, domestic price (lakh Rs per kt), domestic demand share of crude run,
# export discount.
PRODUCTS = [
    ("lpg", 5200.0, 0.02, 0.10),
    ("naphtha", 5000.0, 0.06, 0.04),
    ("gasoline", 6900.0, 0.15, 0.08),
    ("atf", 6700.0, 0.10, 0.06),
    ("diesel", 6600.0, 0.35, 0.06),
    ("fuel_oil", 4100.0, 0.10, 0.05),
    ("bitumen", 3900.0, 0.05, 0.20),
]


def crude_yields(api: float) -> list[float]:
    w = min(1.0, max(0.0, (api - 28.0) / 12.0))
    y = [(1 - w) * h + w * l for h, l in zip(YIELD_HEAVY, YIELD_LIGHT)]
    s = sum(y)
    return [0.99 * v / s for v in y]  # 1% losses


# ------------------------------------------------------------------------------------------
# Planning (LP)


def plan(crudes: int, periods: int, rng: random.Random) -> MpsWriter:
    lp = MpsWriter(f"MRPL_PLAN_{crudes}c_{periods}t", maximize=True)
    slate = CRUDES[:crudes]
    days = 7.0  # Each period is a week.
    cap = {"cdu": 45.0, "vdu": 14.0, "fcc": 6.5, "hcu": 5.8, "reformer": 2.6, "dhds": 11.0}
    op_cost = {"cdu": 60.0, "vdu": 40.0, "fcc": 180.0, "hcu": 260.0, "reformer": 150.0,
               "dhds": 120.0}
    # Conversion yields (fractions of the feed).
    fcc = {"lpg": 0.18, "gasoline": 0.47, "lco": 0.17, "slurry": 0.12}          # 6% coke + gas
    hcu = {"lpg": 0.04, "naphtha": 0.15, "kero": 0.28, "diesel": 0.48}         # 5% gas, H2 gain
    reformer = {"lpg": 0.08, "reformate": 0.84}                                # 8% H2 + gas
    ron = {"light_naphtha": 70.0, "reformate": 98.0, "fcc_gasoline": 92.0}
    gas_sulfur_ppm = {"light_naphtha": 5.0, "reformate": 0.5, "fcc_gasoline": 8.0}
    for t in range(periods):
        demand_scale = rng.uniform(0.9, 1.1)
        for c, (name, api, sulfur, cost) in enumerate(slate):
            price = cost * rng.uniform(0.96, 1.04)
            lp.col(f"buy_{name}_{t}", cost=-price, up=rng.uniform(80, 200))
            lp.col(f"run_{name}_{t}")
            lp.col(f"ctank_{name}_{t}", cost=-2.0, up=120.0)
        for cut in CUTS:
            lp.col(f"sr_{cut}_{t}")
        for u in cap:
            lp.col(f"feed_{u}_{t}", cost=-op_cost[u])
        # Streams.
        for s in ["ln_to_gasoline", "ln_to_naphtha", "hn_to_reformer", "hn_to_naphtha",
                  "vgo_vdu", "vacuum_residue", "light_naphtha", "heavy_naphtha", "gasoil_to_dhds",
                  "kero_to_atf", "kero_to_dhds", "lco_to_dhds", "lco_to_fo", "vgo_to_fcc",
                  "vgo_to_hcu", "residue_to_vdu", "residue_to_fo", "vr_to_bitumen", "vr_to_fo",
                  "slurry_to_fo", "hcu_naphtha", "hcu_kero", "hcu_diesel", "dhds_diesel",
                  "reformate", "fcc_gasoline", "lpg_pool"]:
            lp.col(f"{s}_{t}")
        for name, price, share, discount in PRODUCTS:
            demand = share * cap["cdu"] * days * demand_scale
            lp.col(f"make_{name}_{t}")
            lp.col(f"sell_{name}_{t}", cost=price, up=demand)
            lp.col(f"export_{name}_{t}", cost=price * (1 - discount), up=1.5 * demand)
            lp.col(f"ptank_{name}_{t}", cost=-3.0, up=0.5 * demand)

    for t in range(periods):
        run = {f"run_{n}_{t}": 1.0 for n, *_ in slate}
        for name, api, sulfur, _ in slate:
            terms = {f"buy_{name}_{t}": 1.0, f"run_{name}_{t}": -1.0, f"ctank_{name}_{t}": -1.0}
            if t > 0:
                terms[f"ctank_{name}_{t - 1}"] = 1.0
            lp.row(f"crude_balance_{name}_{t}", terms, 0.0, 0.0)
        lp.row(f"cdu_feed_{t}", {**run, f"feed_cdu_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"cdu_capacity_{t}", {f"feed_cdu_{t}": 1.0}, None, cap["cdu"] * days)
        # High-sulfur metallurgy limit on the CDU blend: average sulfur <= 2.5 %.
        lp.row(f"cdu_sulfur_{t}", {f"run_{n}_{t}": s - 2.5 for n, _, s, _ in slate}, None, 0.0)
        for k, cut in enumerate(CUTS):
            terms = {f"sr_{cut}_{t}": 1.0}
            for name, api, _, _ in slate:
                terms[f"run_{name}_{t}"] = -crude_yields(api)[k]
            lp.row(f"cut_{cut}_{t}", terms, 0.0, 0.0)
        # Gasoil sulfur is far above 10 ppm for every crude: all of it goes through the DHDS.
        lp.row(f"gasoil_route_{t}", {f"sr_gasoil_{t}": 1.0, f"gasoil_to_dhds_{t}": -1.0}, 0.0, 0.0)
        # Naphtha split into light (65%) and heavy (35%) naphtha.
        lp.row(f"naphtha_split_{t}", {f"sr_naphtha_{t}": 1.0, f"light_naphtha_{t}": -1.0,
                                       f"heavy_naphtha_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"naphtha_ratio_{t}", {f"light_naphtha_{t}": 0.65, f"heavy_naphtha_{t}": -0.35},
               0.0, 0.0)
        lp.row(f"kero_split_{t}", {f"sr_kero_{t}": 1.0, f"kero_to_atf_{t}": -1.0,
                                    f"kero_to_dhds_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"residue_split_{t}", {f"sr_residue_{t}": 1.0, f"residue_to_vdu_{t}": -1.0,
                                       f"residue_to_fo_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"vdu_feed_{t}", {f"feed_vdu_{t}": 1.0, f"residue_to_vdu_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"vdu_yield_vgo_{t}", {f"vgo_vdu_{t}": 1.0, f"feed_vdu_{t}": -0.45}, 0.0, 0.0)
        lp.row(f"vdu_yield_vr_{t}", {f"vacuum_residue_{t}": 1.0, f"feed_vdu_{t}": -0.55}, 0.0, 0.0)
        lp.row(f"vgo_split_{t}", {f"sr_vgo_{t}": 1.0, f"vgo_vdu_{t}": 1.0,
                                   f"vgo_to_fcc_{t}": -1.0, f"vgo_to_hcu_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"vr_split_{t}", {f"vacuum_residue_{t}": 1.0, f"vr_to_bitumen_{t}": -1.0,
                                  f"vr_to_fo_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"fcc_feed_{t}", {f"feed_fcc_{t}": 1.0, f"vgo_to_fcc_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"hcu_feed_{t}", {f"feed_hcu_{t}": 1.0, f"vgo_to_hcu_{t}": -1.0}, 0.0, 0.0)
        # Heavy straight-run naphtha and hydrocracker naphtha feed the reformer.
        lp.row(f"hn_split_{t}", {f"heavy_naphtha_{t}": 1.0, f"hn_to_reformer_{t}": -1.0,
                                  f"hn_to_naphtha_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"reformer_feed_{t}", {f"feed_reformer_{t}": 1.0, f"hn_to_reformer_{t}": -1.0,
                                       f"hcu_naphtha_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"fcc_gasoline_yield_{t}", {f"fcc_gasoline_{t}": 1.0,
                                            f"feed_fcc_{t}": -fcc["gasoline"]}, 0.0, 0.0)
        lp.row(f"lco_split_{t}", {f"feed_fcc_{t}": fcc["lco"], f"lco_to_dhds_{t}": -1.0,
                                   f"lco_to_fo_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"slurry_{t}", {f"feed_fcc_{t}": fcc["slurry"], f"slurry_to_fo_{t}": -1.0}, 0.0,
               0.0)
        for s, key in [("hcu_naphtha", "naphtha"), ("hcu_kero", "kero"),
                       ("hcu_diesel", "diesel")]:
            lp.row(f"{s}_yield_{t}", {f"{s}_{t}": 1.0, f"feed_hcu_{t}": -hcu[key]}, 0.0, 0.0)
        lp.row(f"reformate_yield_{t}", {f"reformate_{t}": 1.0,
                                         f"feed_reformer_{t}": -reformer["reformate"]}, 0.0, 0.0)
        # DHDS: gasoil, kerosene and LCO in; 97% out as 8 ppm diesel.
        lp.row(f"dhds_feed_{t}", {f"feed_dhds_{t}": 1.0, f"gasoil_to_dhds_{t}": -1.0,
                                   f"kero_to_dhds_{t}": -1.0, f"lco_to_dhds_{t}": -1.0}, 0.0, 0.0)
        # LCO is hard to treat: at most 20% of the DHDS feed.
        lp.row(f"dhds_lco_limit_{t}", {f"lco_to_dhds_{t}": 1.0, f"feed_dhds_{t}": -0.2}, None, 0.0)
        lp.row(f"dhds_yield_{t}", {f"dhds_diesel_{t}": 1.0, f"feed_dhds_{t}": -0.97}, 0.0, 0.0)
        lp.row(f"lpg_{t}", {f"lpg_pool_{t}": 1.0, f"sr_lpg_{t}": -1.0,
                            f"feed_fcc_{t}": -fcc["lpg"], f"feed_hcu_{t}": -hcu["lpg"],
                            f"feed_reformer_{t}": -reformer["lpg"]}, 0.0, 0.0)
        for u in ["vdu", "fcc", "hcu", "reformer", "dhds"]:
            lp.row(f"{u}_capacity_{t}", {f"feed_{u}_{t}": 1.0}, None, cap[u] * days)

        # Product pools.
        lp.row(f"pool_lpg_{t}", {f"make_lpg_{t}": 1.0, f"lpg_pool_{t}": -1.0}, 0.0, 0.0)
        # Light naphtha: gasoline blendstock or petrochemical naphtha.
        lp.row(f"ln_split_{t}", {f"light_naphtha_{t}": 1.0, f"ln_to_gasoline_{t}": -1.0,
                                  f"ln_to_naphtha_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"pool_naphtha_{t}", {f"make_naphtha_{t}": 1.0, f"ln_to_naphtha_{t}": -1.0,
                                      f"hn_to_naphtha_{t}": -1.0}, 0.0, 0.0)
        gasoline = {"light_naphtha": f"ln_to_gasoline_{t}", "reformate": f"reformate_{t}",
                    "fcc_gasoline": f"fcc_gasoline_{t}"}
        lp.row(f"pool_gasoline_{t}", {f"make_gasoline_{t}": 1.0,
                                       **{v: -1.0 for v in gasoline.values()}}, 0.0, 0.0)
        # BS-VI gasoline: RON >= 91 and sulfur <= 10 ppm (linear blending by mass).
        lp.row(f"spec_gasoline_ron_{t}", {v: ron[k] - 91.0 for k, v in gasoline.items()}, 0.0,
               None)
        lp.row(f"spec_gasoline_sulfur_{t}",
               {v: gas_sulfur_ppm[k] - 10.0 for k, v in gasoline.items()}, None, 0.0)
        lp.row(f"pool_atf_{t}", {f"make_atf_{t}": 1.0, f"kero_to_atf_{t}": -1.0,
                                  f"hcu_kero_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"pool_diesel_{t}", {f"make_diesel_{t}": 1.0, f"dhds_diesel_{t}": -1.0,
                                     f"hcu_diesel_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"pool_fuel_oil_{t}", {f"make_fuel_oil_{t}": 1.0, f"residue_to_fo_{t}": -1.0,
                                       f"vr_to_fo_{t}": -1.0, f"lco_to_fo_{t}": -1.0,
                                       f"slurry_to_fo_{t}": -1.0}, 0.0, 0.0)
        lp.row(f"pool_bitumen_{t}", {f"make_bitumen_{t}": 1.0, f"vr_to_bitumen_{t}": -1.0},
               0.0, 0.0)
        for name, *_ in PRODUCTS:
            terms = {f"make_{name}_{t}": 1.0, f"sell_{name}_{t}": -1.0, f"export_{name}_{t}": -1.0,
                     f"ptank_{name}_{t}": -1.0}
            if t > 0:
                terms[f"ptank_{name}_{t - 1}"] = 1.0
            lp.row(f"product_balance_{name}_{t}", terms, 0.0, 0.0)
    return lp


# ------------------------------------------------------------------------------------------
# Crude receipt scheduling (MILP)


def crude(cargoes: int, days: int, rng: random.Random,
          update: random.Random | None = None) -> MpsWriter:
    lp = MpsWriter(f"MRPL_CRUDE_{cargoes}k_{days}d", maximize=True)
    slate = CRUDES[:6]
    cdu_min, cdu_max = 36.0, 45.0
    tank_cap = {n: 180.0 for n, *_ in slate}
    total_tank = 520.0
    # Netback of each crude (lakh Rs per kt processed): product value minus crude cost; sweet,
    # light crudes are worth more, and the sulfur limit keeps the heavy ones in check.
    netback = {n: 900.0 + 60.0 * (api - 28.0) - 120.0 * s for n, api, s, _ in slate}
    initial = {n: rng.uniform(40, 90) for n, *_ in slate}
    plan_cargo = []
    for k in range(cargoes):
        name, api, sulfur, cost = slate[k % len(slate)]
        size = rng.choice([80.0, 100.0, 130.0])  # Aframax, Suezmax part cargo, Suezmax
        early = rng.randint(0, max(0, days - 4))
        late = min(days - 1, early + rng.randint(2, 4))
        plan_cargo.append((k, name, size, early, late))
    if update is not None:
        # Re-plan the next day: opening stocks as measured (within 5% of the forecast) and one
        # term cargo reporting a delay of one or two days.
        initial = {n: v * update.uniform(0.95, 1.05) for n, v in initial.items()}
        k = update.randrange(cargoes)
        _, name, size, early, late = plan_cargo[k]
        slip = update.randint(1, 2)
        plan_cargo[k] = (k, name, size, min(days - 1, early + slip), min(days - 1, late + slip))
    for k, name, size, early, late in plan_cargo:
        for t in range(early, late + 1):
            # Demurrage after the first two days of the window.
            late_cost = 25.0 * max(0, t - early - 1) * size / 100.0
            lp.col(f"arrive_{k}_{name}_{t}", cost=-late_cost, up=1, integer=True)
        # A term cargo that is not lifted costs a contract penalty (about $6 million).
        lp.col(f"skip_{k}_{name}", cost=-5000.0, up=1, integer=True)
    spot = []
    for s in range(max(2, cargoes // 4)):
        name = slate[rng.randrange(len(slate))][0]
        for t in range(days):
            lp.col(f"spot_{s}_{name}_{t}", cost=-rng.uniform(150, 300), up=1, integer=True)
        spot.append((s, name))
    for t in range(days):
        for n, *_ in slate:
            lp.col(f"tank_{n}_{t}", cost=-0.3, up=tank_cap[n])
            lp.col(f"charge_{n}_{t}", cost=netback[n])
        lp.col(f"cdu_on_{t}", up=1, integer=True)
    for k, name, size, early, late in plan_cargo:
        lp.row(f"lift_{k}", {**{f"arrive_{k}_{name}_{t}": 1.0 for t in range(early, late + 1)},
                             f"skip_{k}_{name}": 1.0}, 1.0, 1.0)
    for t in range(days):
        berth = {}
        for k, name, size, early, late in plan_cargo:
            if early <= t <= late:
                berth[f"arrive_{k}_{name}_{t}"] = 1.0
        for s, name in spot:
            berth[f"spot_{s}_{name}_{t}"] = 1.0
        lp.row(f"berth_{t}", berth, None, 1.0)
        for n, api, sulfur, _ in slate:
            terms = {f"tank_{n}_{t}": -1.0, f"charge_{n}_{t}": -1.0}
            for k, name, size, early, late in plan_cargo:
                if name == n and early <= t <= late:
                    terms[f"arrive_{k}_{name}_{t}"] = size
            for s, name in spot:
                if name == n:
                    terms[f"spot_{s}_{name}_{t}"] = 90.0
            if t > 0:
                terms[f"tank_{n}_{t - 1}"] = 1.0
            rhs = -initial[n] if t == 0 else 0.0
            lp.row(f"tank_balance_{n}_{t}", terms, rhs, rhs)
        charge = {f"charge_{n}_{t}": 1.0 for n, *_ in slate}
        lp.row(f"cdu_max_{t}", {**charge, f"cdu_on_{t}": -cdu_max}, None, 0.0)
        lp.row(f"cdu_min_{t}", {**charge, f"cdu_on_{t}": -cdu_min}, 0.0, None)
        lp.row(f"cdu_sulfur_{t}", {f"charge_{n}_{t}": s - 2.3 for n, _, s, _ in slate}, None, 0.0)
        lp.row(f"tank_farm_{t}", {f"tank_{n}_{t}": 1.0 for n, *_ in slate}, None, total_tank)
        # A shutdown day costs throughput; at most one planned stop in the horizon.
    lp.row("cdu_stops", {f"cdu_on_{t}": 1.0 for t in range(days)}, days - 1, None)
    # Keep a working stock at the end of the horizon.
    lp.row("closing_stock", {f"tank_{n}_{days - 1}": 1.0 for n, *_ in slate}, 150.0, None)
    return lp


# ------------------------------------------------------------------------------------------
# Utility unit commitment (MILP)


def utility(hours: int, rng: random.Random,
            update: random.Random | None = None) -> MpsWriter:
    lp = MpsWriter(f"MRPL_UTILITY_{hours}h")
    # name: max HP steam (t/h), min load fraction, fuel cost per t steam (lakh Rs), no-load cost
    # per hour, startup cost, min up hours.
    boilers = [("boiler1", 150.0, 0.35, 0.021, 0.30, 1.5, 4),
               ("boiler2", 150.0, 0.35, 0.022, 0.30, 1.5, 4),
               ("boiler3", 120.0, 0.40, 0.024, 0.25, 1.2, 4)]
    gtg_mw, gtg_min, gtg_heat_cost, gtg_noload, gtg_start, gtg_up = 35.0, 0.5, 0.075, 0.9, 4.0, 6
    hrsg_steam_per_mw = 2.2  # t/h of HP steam per MW of GT output
    stg_mw_per_t = 0.12      # MW per t/h of HP steam through the back-pressure turbine
    stg_max = 180.0
    tariff = [0.075 if 6 <= h % 24 < 22 else 0.048 for h in range(hours)]  # lakh Rs per MWh
    ramp = 0.3
    for h in range(hours):
        daily = 1.0 + 0.12 * (1 if 8 <= h % 24 < 20 else -1)
        power = 62.0 * daily * rng.uniform(0.97, 1.03)
        hp = 110.0 * daily * rng.uniform(0.95, 1.05)
        mp = 210.0 * daily * rng.uniform(0.95, 1.05)
        if update is not None:
            # Re-plan with the revised demand forecast (each hour within 5% of yesterday's).
            power *= update.uniform(0.95, 1.05)
            hp *= update.uniform(0.95, 1.05)
            mp *= update.uniform(0.95, 1.05)
        for name, cap, mn, fuel, noload, start, _ in boilers:
            lp.col(f"{name}_on_{h}", cost=noload, up=1, integer=True)
            lp.col(f"{name}_start_{h}", cost=start, up=1, integer=True)
            lp.col(f"{name}_steam_{h}", cost=fuel, up=cap)
        lp.col(f"gtg_on_{h}", cost=gtg_noload, up=1, integer=True)
        lp.col(f"gtg_start_{h}", cost=gtg_start, up=1, integer=True)
        lp.col(f"gtg_mw_{h}", cost=gtg_heat_cost, up=gtg_mw)
        lp.col(f"stg_steam_{h}", up=stg_max)
        lp.col(f"letdown_{h}", up=250.0)
        lp.col(f"import_mw_{h}", cost=tariff[h], up=40.0)
        lp.col(f"export_mw_{h}", cost=-0.035, up=10.0)
        lp.col(f"vent_mp_{h}", cost=0.01, up=40.0)
        steam = {f"{n}_steam_{h}": 1.0 for n, *_ in boilers}
        lp.row(f"hp_balance_{h}", {**steam, f"gtg_mw_{h}": hrsg_steam_per_mw,
                                    f"stg_steam_{h}": -1.0, f"letdown_{h}": -1.0}, hp, None)
        lp.row(f"mp_balance_{h}", {f"stg_steam_{h}": 1.0, f"letdown_{h}": 1.0,
                                    f"vent_mp_{h}": -1.0}, mp, None)
        lp.row(f"power_balance_{h}", {f"gtg_mw_{h}": 1.0, f"stg_steam_{h}": stg_mw_per_t,
                                       f"import_mw_{h}": 1.0, f"export_mw_{h}": -1.0}, power,
               None)
        lp.row(f"gtg_max_{h}", {f"gtg_mw_{h}": 1.0, f"gtg_on_{h}": -gtg_mw}, None, 0.0)
        lp.row(f"gtg_min_{h}", {f"gtg_mw_{h}": 1.0, f"gtg_on_{h}": -gtg_min * gtg_mw}, 0.0, None)
        terms = {f"gtg_start_{h}": 1.0, f"gtg_on_{h}": -1.0}
        if h > 0:
            terms[f"gtg_on_{h - 1}"] = 1.0
        lp.row(f"gtg_startup_{h}", terms, 0.0 if h > 0 else -1.0, None)
        for name, cap, mn, *_ in boilers:
            lp.row(f"{name}_max_{h}", {f"{name}_steam_{h}": 1.0, f"{name}_on_{h}": -cap}, None,
                   0.0)
            lp.row(f"{name}_min_{h}", {f"{name}_steam_{h}": 1.0, f"{name}_on_{h}": -mn * cap},
                   0.0, None)
            terms = {f"{name}_start_{h}": 1.0, f"{name}_on_{h}": -1.0}
            if h > 0:
                terms[f"{name}_on_{h - 1}"] = 1.0
            lp.row(f"{name}_startup_{h}", terms, 0.0 if h > 0 else -1.0, None)
            if h > 0:
                lp.row(f"{name}_ramp_up_{h}", {f"{name}_steam_{h}": 1.0,
                                                f"{name}_steam_{h - 1}": -1.0,
                                                f"{name}_start_{h}": -mn * cap},
                       None, ramp * cap)
                lp.row(f"{name}_ramp_down_{h}", {f"{name}_steam_{h - 1}": 1.0,
                                                  f"{name}_steam_{h}": -1.0,
                                                  f"{name}_on_{h}": mn * cap},
                       None, ramp * cap + mn * cap)
    # Minimum up time: a unit started in hour h stays on for the next `up` hours.
    for name, *_, up in boilers + [("gtg", 0, 0, 0, 0, 0, gtg_up)]:
        for h in range(hours):
            for k in range(h + 1, min(hours, h + up)):
                lp.row(f"{name}_minup_{h}_{k}", {f"{name}_on_{k}": 1.0, f"{name}_start_{h}": -1.0},
                       0.0, None)
    return lp


# ------------------------------------------------------------------------------------------

SIZES = {
    "plan": [("small", (4, 4)), ("medium", (6, 13)), ("large", (8, 52))],
    "crude": [("small", (8, 14)), ("medium", (14, 21)), ("large", (22, 30))],
    "utility": [("small", (24,)), ("medium", (48,)), ("large", (96,))],
}
BUILDERS = {"plan": plan, "crude": crude, "utility": utility}


def generate(out: Path, seed: int, update: int | None = None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    for family, sizes in SIZES.items():
        if update is not None and family == "plan":
            continue  # Re-planning variants exist for the scheduling MILPs only.
        for size, args in sizes:
            rng = random.Random(f"{seed}-{family}-{size}")
            extra = [] if update is None else [random.Random(f"{update}-{family}-{size}")]
            model = BUILDERS[family](*args, rng, *extra)
            suffix = "" if update is None else f"_update{update}"
            path = out / f"mrpl_{family}_{size}{suffix}.mps"
            model.write(path)
            n_int = len(model.integer)
            print(f"{path}: {len(model.rows)} rows, {len(model.cols)} columns ({n_int} integer)")


def read_solution(path: Path, *, details: dict | None = None) -> tuple[dict[str, str], dict[str, float], dict[str, float]]:
    """Keep the header/values/duals return; optional details receives all named fields and ranging."""
    header: dict[str, str] = {}
    values: dict[str, float] = {}
    duals: dict[str, float] = {}
    section = None
    fields = []
    data = {"columns": {}, "rows": {}, "ranging": None}
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if not parts:
            continue
        if parts[0] == "#":
            fields = parts[1:]
            continue
        if parts[0].startswith("#"):
            continue
        if parts[0] in ("columns", "rows", "infeasibility_certificate", "unbounded_ray"):
            section = parts[0]
            fields = []
            continue
        if section is None:
            header[parts[0]] = " ".join(parts[1:])
            if parts[0] == "ranging":
                data["ranging"] = {"available": parts[1] == "available", "message": " ".join(parts[2:])}
        elif section in ("columns", "rows"):
            record = dict(zip(fields, parts))
            name = record.pop("name", None)
            if name is None:
                continue
            numbers = {key: float(value) for key, value in record.items()}
            data[section][name] = numbers
            if section == "columns":
                values[name] = numbers.get("value", 0.0)
            elif "dual" in numbers:
                duals[name] = numbers["dual"]
    if details is not None:
        details.update(data)
    return header, values, duals


def report(model: Path, solution: Path) -> None:
    header, x, y = read_solution(solution)
    name = model.stem
    print(f"{name}: status {header.get('status')}, verified {header.get('verified')}, "
          f"objective {float(header.get('objective', 'nan')):,.1f}")

    def series(prefix: str) -> dict[str, list[float]]:
        out: dict[str, list[float]] = {}
        for k, v in x.items():
            if k.startswith(prefix):
                key, t = k[len(prefix):].rsplit("_", 1)
                out.setdefault(key, [])
                idx = int(t)
                while len(out[key]) <= idx:
                    out[key].append(0.0)
                out[key][idx] = v
        return out

    def table(title: str, rows: dict[str, list[float]], fmt: str = "{:7.1f}") -> None:
        rows = {k: v for k, v in rows.items() if any(abs(a) > 1e-6 for a in v)}
        if not rows:
            return
        width = max(len(k) for k in rows)
        periods = max(len(v) for v in rows.values())
        print(f"\n{title}")
        print(" " * (width + 2) + "".join(f"{t + 1:>8}" for t in range(periods)))
        for k, v in sorted(rows.items()):
            print(f"  {k:<{width}}" + "".join(" " + fmt.format(a + 0.0) for a in v))

    if "PLAN" in name.upper():
        table("Crude processed (kt per week)", series("run_"))
        table("Unit feed (kt per week)", series("feed_"))
        table("Domestic sales (kt per week)", series("sell_"))
        table("Exports (kt per week)", series("export_"))
        caps = {k: v for k, v in y.items() if "capacity" in k and abs(v) > 1e-6}
        if caps:
            print("\nMarginal value of capacity (lakh Rs per extra kt), binding constraints only:")
            for k, v in sorted(caps.items()):
                print(f"  {k:<24} {abs(v):9.1f}")
    elif "CRUDE" in name.upper():
        arrivals = sorted((int(k.rsplit("_", 1)[1]), k) for k, v in x.items()
                          if (k.startswith("arrive_") or k.startswith("spot_")) and v > 0.5)
        print("\nBerth schedule:")
        for t, k in arrivals:
            parts = k.split("_")
            kind = "term cargo" if parts[0] == "arrive" else "spot cargo"
            print(f"  day {t + 1:>3}: {kind} {parts[1]} ({'_'.join(parts[2:-1])})")
        skipped = [k for k, v in x.items() if k.startswith("skip_") and v > 0.5]
        if skipped:
            print("  not lifted: " + ", ".join(skipped))
        table("CDU charge by crude (kt per day)", series("charge_"))
        table("Tank stock (kt, end of day)", series("tank_"))
    elif "UTILITY" in name.upper():
        on = series("")
        units = {k[:-3]: v for k, v in on.items() if k.endswith("_on")}
        print("\nUnit commitment (# = on):")
        for unit, v in sorted(units.items()):
            print(f"  {unit:<8} " + "".join("#" if a > 0.5 else "." for a in v))
        table("Steam and power", {k: v for k, v in on.items()
                                  if k.endswith(("_steam", "_mw")) and not k.startswith("stg")},
              "{:6.0f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    g = sub.add_parser("generate")
    g.add_argument("--out", type=Path, default=Path(__file__).parent / "instances")
    g.add_argument("--seed", type=int, default=26119)
    g.add_argument("--update", type=int, default=None,
                   help="also write re-planning variants (same plant, revised data) with this "
                        "seed, named *_update<seed>.mps")
    r = sub.add_parser("report")
    r.add_argument("model", type=Path)
    r.add_argument("solution", type=Path)
    args = parser.parse_args()
    if args.command == "generate":
        generate(args.out, args.seed)
        if args.update is not None:
            generate(args.out, args.seed, args.update)
    else:
        report(args.model, args.solution)


if __name__ == "__main__":
    main()
