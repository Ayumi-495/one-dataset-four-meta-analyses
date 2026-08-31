import { execFileSync } from "node:child_process";
import { existsSync, readFileSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dirname, "..", "..");
const revision = resolve(root, "revision");
const qmdPath = resolve(revision, "tutorial.qmd");
const htmlPath = resolve(revision, "index.html");
const artifactPath = resolve(revision, "Rdata", "study_level_artifacts.rds");
const scriptPath = resolve(revision, "R", "study_level_heterogeneity.R");
const blsmetaProvenancePath = resolve(revision, "data", "blsmeta_sensitivity_verified_outputs.csv");

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
    "blsmeta` is a direct Bayesian study-scale implementation",
    "Historical `blsmeta` sensitivity output",
    "claim-bearing results and are not evidence"
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

function checkInstructional() {
  const qmd = readFileSync(qmdPath, "utf8");
  const required = [
    "## Starting comparator: residual location--scale model",
    "## General multilevel location--scale model",
    "## Categorical random coefficients are a location--scale special case",
    "## Primary full-data categorical model",
    "study-level-metafor-code",
    "study-level-glmmtmb-code",
    "study-level-brms-code",
    "study-level-drmtmb-code",
    "study-level-blsmeta-code",
    "study-level-prior-sensitivity",
    "## Matched-subset direct study-scale bridge",
    "Only six studies contain both fertiliser categories",
    "exp(\\gamma_1/2)",
    "\\operatorname{Var}(u_{jP})",
    "\\mathbf1^\\mathsf{T}\\mathbf V_{\\mathrm{new}}\\mathbf1",
    "predeclared convergence criteria",
    "study-level-blsmeta-historical-results"
  ];
  for (const text of required) assert(qmd.includes(text), `missing instructional material: ${text}`);
  console.log("instructional verification passed");
}

function checkTechnical() {
  const qmd = readFileSync(qmdPath, "utf8");
  const script = readFileSync(scriptPath, "utf8");
  for (const text of [
    "residual_animal <- exp(draws[[draw_column(draws, \"^b_sigma_fertilizeranimal$\")]])",
    "residual_plant <- exp(draws[[draw_column(draws, \"^b_sigma_fertilizerplant$\")]])",
    "median_response_by_fertilizer",
    "boundary = c(boundary, FALSE)"
  ]) assert(script.includes(text), `missing technical correction: ${text}`);
  assert(qmd.includes("separate median sampling variance\nfor each response-by-fertiliser category"),
    "tutorial does not describe response-by-fertiliser sampling-variance scenarios");
  assert(qmd.includes("must be exponentiated before forming SD draws or prediction intervals"),
    "tutorial does not explain the brms log-SD conversion");
  console.log("technical verification passed");
}

function checkPipeline() {
  assert(existsSync(scriptPath), "study-level pipeline script is missing");
  const script = readFileSync(scriptPath, "utf8");
  for (const text of ["measure = \"CVR\"", "correct = TRUE", "fit_brms_category_model", "fit_metafor_category_model", "fit_glmmtmb_category_model", "fit_drmtmb_direct_model", "run_blsmeta_sensitivity", "BLSMETA SENSITIVITY COMPLETED", "saveRDS"]) {
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

function checkBlsmeta() {
  assert(existsSync(blsmetaProvenancePath), "blsmeta provenance table is missing");
  const qmd = readFileSync(qmdPath, "utf8");
  const provenance = readFileSync(blsmetaProvenancePath, "utf8");
  for (const text of [
    "Historical `blsmeta` sensitivity output",
    "Both `blsmeta` fits **failed** the predeclared convergence criteria",
    "claim-bearing results and are not evidence",
    "study-level-blsmeta-historical-results",
    "blsmeta_sensitivity_verified_outputs.csv"
  ]) assert(qmd.includes(text), `missing blsmeta documentation: ${text}`);
  for (const text of [
    "1.11347113985699",
    "2.50497961040809",
    "2.11090635552233",
    "0.854092881471793",
    "failed_convergence",
    "1b759ac",
    "Study-ratio target diagnostics did not fully pass."
  ]) assert(provenance.includes(text), `missing blsmeta provenance value: ${text}`);
  console.log("blsmeta verification passed");
}

function checkRendered() {
  assert(existsSync(htmlPath), "revision/index.html is missing");
  const html = readFileSync(htmlPath, "utf8");
  for (const text of ["Where does heterogeneity occur?", "General multilevel location", "Categorical random coefficients", "blsmeta", "Historical <code>blsmeta</code> sensitivity output", "FAILED: study-ratio target criteria not met", "1.11 [0.89, 1.40]", "2.50 [1.00, 6.60]", "2.11 [1.47, 3.03]", "0.85 [0.20, 3.84]", "Posterior plant-to-animal SD ratios", "2.18", "1.67"]) {
    assert(html.includes(text), `rendered tutorial is missing ${text}`);
  }
  for (const forbidden of ["lower stability", "strongly governs the predictability", "stability gap"]) {
    assert(!html.includes(forbidden), `stale rendered interpretation remains: ${forbidden}`);
  }
  console.log("rendered verification passed");
}

const checks = { scope: checkScope, source: checkSource, instructional: checkInstructional, technical: checkTechnical, pipeline: checkPipeline, blsmeta: checkBlsmeta, rendered: checkRendered };
const selected = process.argv[2];
if (!checks[selected]) fail(`choose one of: ${Object.keys(checks).join(", ")}`);
checks[selected]();
