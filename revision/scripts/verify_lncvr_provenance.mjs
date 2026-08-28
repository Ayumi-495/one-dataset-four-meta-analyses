#!/usr/bin/env node
import { execFileSync } from "node:child_process";
import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";

const revisionRoot = resolve(import.meta.dirname, "..");
const repositoryRoot = resolve(revisionRoot, "..");
const historicCommit = "1b759ac4cc00c5959b172462ebb2acad5f5ca6d2";
const historicResult = [
  "1.05129579235815",
  "0.0991992778030649",
  "8.63061848467873"
];

function gitShow(path) {
  return execFileSync("git", ["show", `${historicCommit}:${path}`], {
    cwd: repositoryRoot, encoding: "utf8"
  });
}

function requireMatch(text, pattern, label) {
  if (!pattern.test(text)) throw new Error(`Missing ${label}`);
}

function historical() {
  const csv = gitShow("results/study_level_heterogeneity/model_summary.csv");
  const script = gitShow("R/study_level_heterogeneity/full_bayesian.R");
  historicResult.forEach((value) => requireMatch(csv, new RegExp(value.replaceAll(".", "\\.")), value));
  requireMatch(script, /study_sd_prior <- if \(sensitivity\) "exponential\(1\)" else "exponential\(2\)"/,
    "historical study-SD prior selection");
  requireMatch(script, /set_prior\(study_sd_prior, class = "sd", group = "study_ID"\)/,
    "historical study-SD prior application");
  requireMatch(script, /full_iter\s*=\s*2000/, "historical iterations");
  requireMatch(script, /full_warmup\s*=\s*1000/, "historical warmup");
  requireMatch(script, /backend\s*=\s*"cmdstanr"/, "historical backend");
  console.log("historical provenance verification passed");
}

function configuration() {
  const script = readFileSync(resolve(revisionRoot, "R/study_level_heterogeneity.R"), "utf8");
  requireMatch(script, /historical_primary_spec/, "named recovered specification");
  requireMatch(script, /set_prior\("exponential\(2\)", class = "sd", group = "study_ID"\)/,
    "explicit recovered study-SD prior");
  requireMatch(script, /set_prior\("normal\(0, 1\)", class = "b"\)/, "location prior");
  requireMatch(script, /set_prior\("normal\(-1, 1\)", class = "b", dpar = "sigma"\)/,
    "residual-scale prior");
  requireMatch(script, /backend\s*=\s*specification\$backend/, "explicit backend");
  requireMatch(script, /iter\s*=\s*specification\$iter/, "explicit iterations");
  requireMatch(script, /correct\s*=\s*TRUE/, "explicit CVR correction");
  console.log("configuration verification passed");
}

function reproduction() {
  const artifact = resolve(revisionRoot, "Rdata/study_level_artifacts.rds");
  if (!existsSync(artifact)) throw new Error(`Missing fitted artifact: ${artifact}`);
  const expression = String.raw`
    a <- readRDS("${artifact.replaceAll("\\", "\\\\").replaceAll('"', '\\"')}")
    stopifnot(identical(a$model_provenance$primary_specification, "recovered_commit_1b759ac"))
    fit <- a$brms_fits$lnCVR
    p <- brms::prior_summary(fit)
    stopifnot(any(p$prior == "exponential(2)" & p$class == "sd" & p$group == "study_ID"))
    x <- subset(a$ratios, response == "lnCVR" & component == "study_sd_ratio")
    stopifnot(nrow(x) == 1L, all(is.finite(unlist(x[c("estimate", "lower", "upper")]))) )
    old <- c(1.05129579235815, 0.0991992778030649, 8.63061848467873)
    now <- unlist(x[c("estimate", "lower", "upper")])
    cat(sprintf("historical reproduction verification passed; reproduced=%.12f [%.12f, %.12f]; delta=%.12f [%.12f, %.12f]\\n", now[1], now[2], now[3], now[1]-old[1], now[2]-old[2], now[3]-old[3]))
  `;
  process.stdout.write(execFileSync("Rscript", ["-e", expression], { cwd: revisionRoot, encoding: "utf8" }));
}

const action = process.argv[2];
if (action === "historical") historical();
else if (action === "configuration") configuration();
else if (action === "reproduction") reproduction();
else throw new Error("Usage: verify_lncvr_provenance.mjs <historical|configuration|reproduction>");
