"use client";

import { useState } from "react";
import { motion, AnimatePresence } from "framer-motion";
import { FlaskConical, Atom, Activity, ArrowRight, RefreshCcw, Info, Dna } from "lucide-react";
import { clsx, type ClassValue } from "clsx";
import { twMerge } from "tailwind-merge";

function cn(...inputs: ClassValue[]) {
  return twMerge(clsx(inputs));
}

const DEFAULT_DRUGS = {
  drug1: "CCO", // Ethanol
  drug2: "CN",  // Methylamine
};

interface DrugResult {
  smiles: string;
  iupac: string;
  description?: string;
  image: string;
}

interface PredictionResult {
  drug1: DrugResult;
  drug2: DrugResult;
  prediction: {
    label: string;
    probability: number;
  };
}

export default function Home() {
  const [smiles1, setSmiles1] = useState("");
  const [smiles2, setSmiles2] = useState("");
  const [loading, setLoading] = useState(false);
  const [result, setResult] = useState<PredictionResult | null>(null);
  const [error, setError] = useState("");

  const handlePredict = async () => {
    if (!smiles1.trim() || !smiles2.trim()) {
      setError("Please enter both SMILES strings.");
      return;
    }
    setError("");
    setLoading(true);
    setResult(null);

    try {
      const response = await fetch("http://localhost:8000/predict", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ smiles1, smiles2 }),
      });

      if (!response.ok) throw new Error("Failed to fetch interaction data");

      const data = await response.json();
      if (data.error) throw new Error(data.error);

      setResult(data);
    } catch (err: any) {
      setError(err.message || "Something went wrong.");
    } finally {
      setLoading(false);
    }
  };

  const loadDefaults = () => {
    setSmiles1(DEFAULT_DRUGS.drug1);
    setSmiles2(DEFAULT_DRUGS.drug2);
    setError("");
  };

  return (
    <div className="min-h-screen bg-[#0a0a0a] text-zinc-100 selection:bg-cyan-500/30 selection:text-cyan-200 overflow-x-hidden font-sans">
      <div className="fixed inset-0 -z-10 h-full w-full bg-[radial-gradient(ellipse_80%_80%_at_50%_-20%,rgba(120,119,198,0.15),rgba(255,255,255,0))]"></div>

      <main className="container mx-auto px-4 py-16 md:py-24 max-w-5xl">
        <motion.div
          initial={{ opacity: 0, y: 20 }}
          animate={{ opacity: 1, y: 0 }}
          className="text-center mb-16 space-y-4"
        >
          <div className="inline-flex items-center justify-center p-2 bg-zinc-900/50 rounded-full border border-zinc-800 mb-4 backdrop-blur-sm">
            <span className="flex h-2 w-2 rounded-full bg-cyan-500 mr-2 animate-pulse"></span>
            <span className="text-xs font-medium text-cyan-500 tracking-wider uppercase">GNN Powered Analysis</span>
          </div>
          <h1 className="text-4xl md:text-6xl font-bold tracking-tight bg-clip-text text-transparent bg-gradient-to-b from-white to-zinc-400">
            Drug Interaction Predictor
          </h1>
          <p className="text-zinc-400 max-w-2xl mx-auto text-lg">
            Analyze molecular compatibility and predict potential interactions using advanced Graph Neural Networks.
          </p>
        </motion.div>

        <div className="grid lg:grid-cols-12 gap-8">
          {/* Input Section */}
          <motion.div
            initial={{ opacity: 0, x: -20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.2 }}
            className="lg:col-span-4 space-y-6"
          >
            <div className="p-6 rounded-2xl bg-zinc-900/40 border border-zinc-800/50 backdrop-blur-xl shadow-2xl">
              <div className="flex items-center justify-between mb-6">
                <h2 className="text-xl font-semibold flex items-center gap-2">
                  <FlaskConical className="w-5 h-5 text-cyan-400" />
                  Configurations
                </h2>
                <button
                  onClick={loadDefaults}
                  className="text-xs text-zinc-500 hover:text-cyan-400 transition-colors flex items-center gap-1"
                >
                  <RefreshCcw className="w-3 h-3" /> Auto-fill
                </button>
              </div>

              <div className="space-y-4">
                <div className="group">
                  <label className="text-xs font-medium text-zinc-500 uppercase tracking-wider mb-1.5 block group-focus-within:text-cyan-500 transition-colors">Molecule A (SMILES)</label>
                  <input
                    type="text"
                    value={smiles1}
                    onChange={(e) => setSmiles1(e.target.value)}
                    placeholder="e.g. CCO"
                    className="w-full bg-black/40 border border-zinc-800 rounded-lg px-4 py-3 text-sm focus:outline-none focus:ring-2 focus:ring-cyan-500/20 focus:border-cyan-500/50 transition-all font-mono"
                  />
                </div>

                <div className="flex justify-center -my-2 relative z-10">
                  <div className="bg-zinc-900/80 p-1.5 rounded-full border border-zinc-800 text-zinc-500">
                    <ArrowRight className="w-4 h-4 rotate-90" />
                  </div>
                </div>

                <div className="group">
                  <label className="text-xs font-medium text-zinc-500 uppercase tracking-wider mb-1.5 block group-focus-within:text-purple-500 transition-colors">Molecule B (SMILES)</label>
                  <input
                    type="text"
                    value={smiles2}
                    onChange={(e) => setSmiles2(e.target.value)}
                    placeholder="e.g. CN"
                    className="w-full bg-black/40 border border-zinc-800 rounded-lg px-4 py-3 text-sm focus:outline-none focus:ring-2 focus:ring-purple-500/20 focus:border-purple-500/50 transition-all font-mono"
                  />
                </div>
              </div>

              {error && (
                <div className="mt-4 p-3 rounded-lg bg-red-500/10 border border-red-500/20 text-red-400 text-sm flex items-start gap-2">
                  <Info className="w-4 h-4 mt-0.5 shrink-0" />
                  {error}
                </div>
              )}

              <button
                onClick={handlePredict}
                disabled={loading}
                className="w-full mt-6 bg-gradient-to-r from-cyan-600 to-blue-600 hover:from-cyan-500 hover:to-blue-500 text-white font-medium py-3 rounded-lg transition-all transform active:scale-[0.98] disabled:opacity-50 disabled:cursor-not-allowed shadow-lg shadow-cyan-900/20 flex items-center justify-center gap-2"
              >
                {loading ? (
                  <>
                    <RefreshCcw className="w-4 h-4 animate-spin" /> Analyzing...
                  </>
                ) : (
                  <>
                    <Activity className="w-4 h-4" /> Analyze Interaction
                  </>
                )}
              </button>
            </div>
          </motion.div>

          {/* Result Section */}
          <motion.div
            initial={{ opacity: 0, x: 20 }}
            animate={{ opacity: 1, x: 0 }}
            transition={{ delay: 0.3 }}
            className="lg:col-span-8"
          >
            <div className="h-full rounded-2xl bg-zinc-900/20 border border-zinc-800/30 backdrop-blur-sm p-1">
              {result ? (
                <div className="h-full flex flex-col">
                  {/* Verdict Header */}
                  <div className={cn(
                    "p-6 rounded-t-xl border-b border-zinc-800/50 bg-gradient-to-r",
                    result.prediction.probability > 0.5
                      ? "from-red-900/20 to-orange-900/10 border-red-900/30"
                      : "from-green-900/20 to-emerald-900/10 border-green-900/30"
                  )}>
                    <div className="flex items-center justify-between">
                      <div>
                        <p className="text-zinc-400 text-sm font-medium uppercase tracking-wider mb-1">Analysis Result</p>
                        <h3 className={cn("text-3xl font-bold flex items-center gap-3",
                          result.prediction.probability > 0.5 ? "text-red-400" : "text-green-400"
                        )}>
                          {result.prediction.probability > 0.5 ? (
                            <>High Interaction Probability</>
                          ) : (
                            <>Low Interaction Probability</>
                          )}
                        </h3>
                      </div>
                      <div className="text-right">
                        <span className="text-4xl font-mono font-bold tracking-tighter opacity-80">
                          {(result.prediction.probability * 100).toFixed(1)}%
                        </span>
                        <p className="text-xs text-zinc-500 mt-1">Confidence Score</p>
                      </div>
                    </div>
                  </div>

                  {/* Molecules Display */}
                  <div className="flex-1 p-6 grid md:grid-cols-2 gap-6">
                    <MoleculeCard data={result.drug1} label="Molecule A" color="cyan" />
                    <MoleculeCard data={result.drug2} label="Molecule B" color="purple" />
                  </div>
                </div>
              ) : (
                <div className="h-full min-h-[400px] flex flex-col items-center justify-center text-zinc-600 space-y-4">
                  <Dna className="w-16 h-16 opacity-20" />
                  <p className="text-lg font-medium">Ready to analyze</p>
                  <p className="text-sm max-w-sm text-center opacity-60">Enter the SMILES strings for two molecules to generate structural visualizations and interaction predictions.</p>
                </div>
              )}
            </div>
          </motion.div>
        </div>
      </main>
    </div>
  );
}

function MoleculeCard({ data, label, color }: { data: DrugResult, label: string, color: "cyan" | "purple" }) {
  return (
    <div className="bg-black/20 rounded-xl border border-zinc-800/50 overflow-hidden flex flex-col">
      <div className={cn("px-4 py-3 border-b border-zinc-800/50 backdrop-blur-sm",
        color === "cyan" ? "bg-cyan-950/10" : "bg-purple-950/10"
      )}>
        <span className={cn("text-xs font-bold uppercase tracking-wider",
          color === "cyan" ? "text-cyan-400" : "text-purple-400"
        )}>{label}</span>
      </div>

      <div className="flex-1 flex items-center justify-center p-6 bg-white/5 relative group min-h-[200px]">
        {/* Grid pattern background */}
        <div className="absolute inset-0 bg-[linear-gradient(to_right,#80808012_1px,transparent_1px),linear-gradient(to_bottom,#80808012_1px,transparent_1px)] bg-[size:14px_14px]"></div>

        {data.image ? (
          <img
            src={`data:image/png;base64,${data.image}`}
            alt={label}
            className="relative z-10 max-h-48 object-contain mix-blend-screen drop-shadow-2xl transition-transform group-hover:scale-105 duration-500"
          />
        ) : (
          <div className="relative z-10 w-full h-40 flex items-center justify-center text-zinc-700">
            <Atom className="w-12 h-12 opacity-50" />
          </div>
        )}
      </div>

      <div className="p-4 space-y-4 bg-zinc-900/30">
        <div>
          <label className="text-[10px] uppercase text-zinc-500 font-semibold tracking-wider">IUPAC Name</label>
          <p className="text-sm text-zinc-300 font-medium leading-snug break-words">{data.iupac || "Unknown"}</p>
        </div>

        <div>
          <label className="text-[10px] uppercase text-zinc-500 font-semibold tracking-wider flex items-center gap-1.5 mb-1.5">
            <Info className="w-3 h-3" /> Description & Context
          </label>
          <div className="max-h-48 overflow-y-auto pr-2 scrollbar-thin scrollbar-thumb-zinc-700 scrollbar-track-transparent">
            <p className="text-sm text-zinc-300 leading-relaxed">
              {data.description || "No description available for this molecule."}
            </p>
          </div>
        </div>

        <div>
          <label className="text-[10px] uppercase text-zinc-500 font-semibold tracking-wider">SMILES</label>
          <p className="text-[10px] text-zinc-600 font-mono break-all line-clamp-1" title={data.smiles}>{data.smiles}</p>
        </div>
      </div>
    </div>
  );
}
