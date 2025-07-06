/* frontend/app.js  –  clean, single-responsibility */

const form   = document.getElementById("form");
const smiles = document.getElementById("smiles");
const out    = document.getElementById("out");

const result = document.getElementById("result");   // ↘ table + charts
const tbl    = document.getElementById("tbl");
const gauge  = document.getElementById("gauge");
const radar  = document.getElementById("radar");

/* ---------- helpers -------------------------------------------------- */
function resetUI() {
  result.style.display = "none";
  gauge.innerHTML = radar.innerHTML = tbl.innerHTML = "";
}

function showMessage(msg, isErr = false) {
  resetUI();
  out.className = isErr ? "error" : "";
  out.textContent = msg;
}

function buildTable(scores) {
  tbl.innerHTML = "<tr><th>Parameter</th><th>Score (0-1)</th></tr>";
  Object.entries(scores).forEach(([k, v]) => {
    tbl.insertAdjacentHTML(
      "beforeend",
      `<tr><td>${k}</td><td>${v.toFixed(3)}</td></tr>`
    );
  });
}

function buildGauge(eri) {
  const data = [
    {
      type: "indicator",
      mode: "gauge+number",
      value: eri,
      number: { suffix: " / 100" },
      gauge: {
        axis: { range: [0, 100] },
        bar: { color: "#10b981" },
        steps: [
          { range: [0, 33],  color: "#d1fae5" },
          { range: [33, 66], color: "#fef9c3" },
          { range: [66, 100], color: "#fee2e2" }
        ]
      }
    }
  ];
  Plotly.newPlot(gauge, data, { margin: { t: 0, b: 0, l: 0, r: 0 } }, { displayModeBar: false });
}

function buildRadar(scores) {
  const labels = Object.keys(scores);
  const vals   = Object.values(scores).map(v => +v.toFixed(3));
  Plotly.newPlot(
    radar,
    [{
      type: "scatterpolar",
      r:    [...vals, vals[0]],
      theta:[...labels, labels[0]],
      fill: "toself"
    }],
    { polar: { radialaxis: { range: [0, 1] } }, margin: { t: 30 } },
    { displayModeBar: false }
  );
}

/* ---------- main handler --------------------------------------------- */
form.addEventListener("submit", async (e) => {
  e.preventDefault();

  const s = smiles.value.trim();
  if (!s) { showMessage("⚠️ Enter a SMILES string.", true); return; }

  showMessage("⏳ Calculating…");

  try {
    const r = await fetch("/api/assess", {
      method:  "POST",
      headers: { "Content-Type": "application/json" },
      body:    JSON.stringify({ smiles: s })
    });
    if (!r.ok) throw new Error((await r.json()).detail || r.statusText);

    const data = await r.json();

    /* build UI */
    buildTable(data.scores);
    buildGauge(data.eri);
    buildRadar(data.scores);

    out.textContent = "";
    result.style.display = "block";
  } catch (err) {
    showMessage("❌ " + err.message, true);
  }
});
