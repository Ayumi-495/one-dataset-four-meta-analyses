import { execFileSync } from "node:child_process";
import { existsSync, readFileSync, statSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dirname, "..", "..");
const revision = resolve(root, "revision");
const qmdPath = resolve(revision, "tutorial.qmd");
const htmlPath = resolve(revision, "index.html");
const artifactPath = resolve(revision, "Rdata", "study_level_artifacts.rds");
const scriptPath = resolve(revision, "R", "study_level_heterogeneity.R");

function fail(message) {
  console.error(`VERIFY FAILED: ${message}`);
  process.exit(1);
}

function command(command, args) {
  return execFileSync(command, args, { cwd: root, encoding: "utf8" });
}

function assert(condition, message) {
  if (!condition) fail(message);
}

function checkScope() {
  const changed = command("git", ["diff", "--name-only"])
    .split("\n").filter(Boolean)
    .concat(command("git", ["ls-files", "--others", "--exclude-standard"])
      .split("\n").filter(Boolean));
  const outOfScope = changed.filter((path) => path !== "GATES.md" && !path.startsWith("revision/"));
  assert(outOfScope.length === 0, `out-of-scope working-tree paths: ${outOfScope.join(", ")}`);
  command("git", ["diff", "--quiet", "--", "tutorial.qmd", "index.html"]);
  console.log("scope verification passed");
}

function checkSource() {
  const qmd = readFileSync(qmdPath, "utf8");
  const required = [
    "# Extension: Where does heterogeneity occur?",
    "correct = TRUE",
    "Standardised mean differences, odds ratios, and Fisher's z-transformed correlations",
    "does not directly estimate a between-group variability contrast",
    "Full categorical model",
    "Direct `drmTMB` study-scale model",
    "`metafor` bridge",
    "observed new comparison",
    "mean of $k$ new comparisons",
    "does not by itself establish poorer transferability",
    "blsmeta` is an optional Bayesian comparator"
  ];
  for (const text of required) assert(qmd.includes(text), `missing required source text: ${text}`);
  const cvrStarts = [...qmd.matchAll(/escalc\(measure\s*=\s*"CVR"/g)].map((match) => match.index);
  assert(cvrStarts.length > 0 && cvrStarts.every((start) => qmd.slice(start, start + 700).includes("correct = TRUE")),
    "each visible CVR calculation must explicitly use correct = TRUE");
  for (const forbidden of ["lower stability", "strongly governs the predictability", "stability gap", "evidence is far less transferable"]) {
    assert(!qmd.includes(forbidden), `stale interpretation remains: ${forbidden}`);
  }
  assert(!qmd.includes("read.csv(text ="), "embedded stale result tables remain in source");
  console.log("source verification passed");
}

function checkPipeline() {
  assert(existsSync(scriptPath), "study-level pipeline script is missing");
  const script = readFileSync(scriptPath, "utf8");
  for (const text of ["measure = \"CVR\"", "correct = TRUE", "fit_brms_category_model", "fit_metafor_category_model", "saveRDS"]) {
    assert(script.includes(text), `pipeline does not contain ${text}`);
  }
  assert(existsSync(artifactPath), "regenerated study-level artifact is missing");
  const check = [
    "x<-readRDS('Rdata/study_level_artifacts.rds')",
    "stopifnot(identical(unname(unlist(x$counts)),c(318L,36L,6L,86L,232L,30L)))",
    "stopifnot(all(c('ratios','predictions','diagnostics','frequentist','matched_data') %in% names(x)))",
    "stopifnot(all(x$diagnostics$max_rhat <= 1.01), all(x$diagnostics$divergences == 0))",
    "cat('artifact verification passed\\n')"
  ].join(";");
  execFileSync("Rscript", ["-e", check], { cwd: revision, stdio: "inherit" });
  console.log("pipeline verification passed");
}

function checkRendered() {
  assert(existsSync(htmlPath), "revision/index.html is missing");
  assert(statSync(htmlPath).mtimeMs >= statSync(qmdPath).mtimeMs,
    "revision/index.html is older than revision/tutorial.qmd");
  const html = readFileSync(htmlPath, "utf8");
  for (const text of ["Where does heterogeneity occur?", "Posterior plant-to-animal SD ratios", "2.18", "1.15", "1.67"]) {
    assert(html.includes(text), `rendered tutorial is missing ${text}`);
  }
  for (const forbidden of ["lower stability", "strongly governs the predictability", "stability gap"]) {
    assert(!html.includes(forbidden), `stale rendered interpretation remains: ${forbidden}`);
  }
  console.log("rendered verification passed");
}

const checks = { scope: checkScope, source: checkSource, pipeline: checkPipeline, rendered: checkRendered };
const selected = process.argv[2];
if (!checks[selected]) fail(`choose one of: ${Object.keys(checks).join(", ")}`);
checks[selected]();
