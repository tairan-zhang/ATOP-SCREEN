import {Button} from "@atop/ui";
import { CenterMorphModal, CenterMorphModalContent } from "@atop/ui/motion/center-morph-modal";
import { useState } from "react";
import type { Project } from "../types/atop";
import { fileName } from "../domain/format";
import { Icon } from "./Icon";

function formatDate(value: string) {
  const date = new Date(value);
  if (Number.isNaN(date.getTime())) return "—";
  return new Intl.DateTimeFormat(undefined, { dateStyle: "medium", timeStyle: "short" }).format(date);
}

export function ProjectLibrary({ busyId, onClose, onDelete, onOpen, onOpenAnalysis, onOpenGsea, projects }: {
  busyId: string;
  onClose: () => void;
  onDelete: (project: Project) => void;
  onOpen: (project: Project) => void;
  onOpenAnalysis: (project: Project, analysisId: string) => void;
  onOpenGsea: (project: Project, analysisId: string, gseaId: string) => void;
  projects: Project[];
}) {
  const [confirmDeleteId, setConfirmDeleteId] = useState("");
  const [expandedProjectIds, setExpandedProjectIds] = useState<Set<string>>(() => new Set());

  function toggleProject(projectId: string) {
    setExpandedProjectIds((current) => {
      const next = new Set(current);
      if (next.has(projectId)) next.delete(projectId);
      else next.add(projectId);
      return next;
    });
  }

  return (
    <CenterMorphModal open onOpenChange={next => {if (!next) onClose();}}><CenterMorphModalContent ariaLabel="Manage projects" className="screen-modal">
      <section className="project-library">
        <header className="project-library__header">
          <div><h2 id="project-library-title">Manage projects</h2></div>
        </header>
        <div className="project-library__list">
          {projects.length ? projects.map((project) => {
            const analyses = project.analyses ?? [];
            const expanded = expandedProjectIds.has(project.id);
            return <section className={`project-library__group${expanded ? " is-expanded" : ""}`} key={project.id}>
              <article className="project-library__row">
                <span className="project-library__icon"><Icon name="folder" /></span>
                <div className="project-library__main">
                  
                  <div><strong>{project.name}</strong><span className={`project-status project-status--${project.status}`}>{project.status}</span></div>
                  <span>{project.inputFileName || "No screen data yet"} · {formatDate(project.updatedAt)}</span>
                  <small title={project.workspacePath}>{fileName(project.workspacePath)}</small>
                  <Button variant="secondary" aria-controls={`project-analyses-${project.id}`} aria-expanded={expanded} className="project-library__disclosure" onClick={() => toggleProject(project.id)} type="button"><Icon name="chevronDown" size={14} /><span>{analyses.length} saved Screen {analyses.length === 1 ? "analysis" : "analyses"}</span></Button>
                </div>
                <div className="project-library__actions">
                  <Button className="primary-button" disabled={busyId === project.id} onClick={() => onOpen(project)} type="button"><Icon name="folder" />Open project</Button>
                  {confirmDeleteId === project.id ? (
                    <div className="project-delete-confirm"><span>Remove from list? Folder files stay.</span><Button variant="secondary" disabled={busyId === project.id} onClick={() => onDelete(project)} type="button">Remove</Button><Button variant="secondary" onClick={() => setConfirmDeleteId("")} type="button">Cancel</Button></div>
                  ) : <Button variant="secondary" aria-label={`Remove ${project.name} from the project list`} className="project-delete-button" disabled={busyId === project.id} onClick={() => setConfirmDeleteId(project.id)} type="button"><Icon name="trash" size={15} /></Button>}
                </div>
              </article>
              {expanded ? <div className="project-library__analyses" id={`project-analyses-${project.id}`}>
                {analyses.length ? analyses.map((analysis) => <article className="project-library__analysis" key={analysis.id}>
                  <div className="project-library__analysis-main"><span className="project-library__branch-icon"><Icon name="flask" size={15} /></span><div><strong>{analysis.name}</strong><small>{analysis.geneRows.toLocaleString()} gene rows · {formatDate(analysis.completedAt ?? project.updatedAt)}</small></div></div>
                  <div className="project-library__analysis-actions"><Button variant="secondary" className="secondary-button" disabled={busyId === project.id} onClick={() => onOpenAnalysis(project, analysis.id)} type="button">Open results</Button></div>
                  {analysis.gseaRuns.length > 0 && <div className="project-library__gsea-runs">
                    <div className="project-library__gsea-label">GSEA runs</div>
                    {analysis.gseaRuns.map((run) => <Button variant="secondary" className="project-library__gsea-run" disabled={busyId === project.id} key={run.id} onClick={() => onOpenGsea(project, analysis.id, run.id)} type="button"><Icon name="bar" size={14} /><span><strong>{run.name}</strong><small>{run.pathwayRows.toLocaleString()} pathway rows · {formatDate(run.completedAt ?? project.updatedAt)}</small></span><Icon name="chevronRight" size={14} /></Button>)}
                  </div>}
                </article>) : <div className="project-library__analyses-empty">No completed Screen analyses in this project.</div>}
              </div> : null}
            </section>;
          }) : <div className="project-library__empty"><Icon name="folder" size={25} /><strong>No saved projects yet</strong></div>}
        </div>
      </section>
    </CenterMorphModalContent></CenterMorphModal>
  );
}
