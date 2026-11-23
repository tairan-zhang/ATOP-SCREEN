import type { View } from "../../types/atop";
import { QuickAction } from "../../components/QuickAction";

type HomeProps = { projectName: string; analysisReady: boolean; onNavigate: (view: View) => void };

export function HomeView({ projectName, analysisReady, onNavigate }: HomeProps) {
  return (
    <div className="home-view">
      <h2 className="project-heading">{projectName}</h2>
      <div className="quick-action-grid">
        <QuickAction icon="flask" label="Analyze counts" onClick={() => onNavigate("analysis")} />
        <QuickAction icon="bar" label="Enrich pathways" onClick={() => onNavigate("gsea")}
          disabled={!analysisReady} disabledReason="Run screen analysis first" />
        <QuickAction icon="image" label="Plot results" onClick={() => onNavigate("visualize")}
          disabled={!analysisReady} disabledReason="Run screen analysis first" />
      </div>
    </div>
  );
}
