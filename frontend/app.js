// frontend/app.js
document.getElementById("form").addEventListener("submit", async (e) => {
    e.preventDefault();
    const smiles = document.getElementById("smiles").value.trim();
    const r = await fetch("/api/assess", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ smiles }),
    });
    const data = await r.json();
    document.getElementById("out").textContent = JSON.stringify(data, null, 2);
  });
