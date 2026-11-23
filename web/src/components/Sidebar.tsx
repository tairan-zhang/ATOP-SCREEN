import { useRef } from "react";
import { Button } from "@atop/ui";
import { AnimatedSidebar, AnimatedSidebarHeader, AnimatedSidebarContent, AnimatedSidebarGroup, AnimatedSidebarGroupLabel } from "@atop/ui/motion/animated-sidebar";
import type { ProcessResult, GseaResult } from "@atop-screen/ipc-contracts";
import type { PlotSection, Project, View } from "../types/atop";
import { WorkspaceResizeHandle } from "../workspace-resize";
import { Icon } from "./Icon";

type ModuleLink = { id: string; label: string; icon: string; view: View; plotSection?: PlotSection; requiresAnalysis?: boolean; requiresGsea?: boolean };

const modules: { id: string; label: string; icon: string; children: ModuleLink[] }[] = [
  { id: "screen", label: "Screen analysis", icon: "flask", children: [
    { id: "screen-setup", label: "Analysis setup", icon: "activity", view: "analysis" },
    { id: "screen-results", label: "Results", icon: "database", view: "results", requiresAnalysis: true },
  ] },
  { id: "pathway", label: "Pathway enrichment", icon: "bar", children: [
    { id: "pathway-setup", label: "Analysis setup", icon: "activity", view: "gsea", requiresAnalysis: true },
    { id: "pathway-results", label: "Results", icon: "database", view: "gseaResults", requiresGsea: true },
  ] },
  { id: "visualization", label: "Visualization", icon: "image", children: [
    { id: "visual-sgrna", label: "sgRNA", icon: "activity", view: "visualize", plotSection: "sgrna", requiresAnalysis: true },
    { id: "visual-gene", label: "Gene", icon: "image", view: "visualize", plotSection: "gene", requiresAnalysis: true },
    { id: "visual-pathway", label: "Pathway", icon: "bar", view: "visualize", plotSection: "pathway", requiresAnalysis: true },
  ] },
];

type Props = {
  project: Project; view: View; plotSection: PlotSection; analysis: ProcessResult | null; gsea: GseaResult | null;
  onCreateProject: () => void; onNewAnalysis: () => void; onOpenLibrary: () => void;
  onNavigate: (view: View) => void; onPlotSection: (section: PlotSection) => void;
  onOpenAnalysis: (id: string) => void; onOpenGsea: (id: string, runId: string) => void;
  onWidth: (width: number) => void; open: boolean;
};

export function Sidebar(props: Props) {
  const sidebar = useRef<HTMLElement>(null);
  return <AnimatedSidebar ref={sidebar} id="workspace-project-sidebar" collapsible="offcanvas"
    ariaLabel="Project navigation" panelClassName="workspace-sidebar" inert={!props.open}>
    <AnimatedSidebarHeader className="workspace-sidebar-header">
      <div className="sidebar-brand"><div className="brand-mode-switch"><img src={import.meta.env.BASE_URL + "brand/atop-glass.png"} alt=""/> ATOP <span className="muted">/ SCREEN</span></div></div>
      <Button className="product-page-nav" variant="ghost" aria-pressed={props.view === "home"} onClick={() => props.onNavigate("home")}><Icon name="home" />Home</Button>
      <Button className="new-session-button" variant="ghost" onClick={props.onNewAnalysis}><Icon name="newAnalysis"/>New analysis</Button>
      <div className="sidebar-toolbar"><Button variant="ghost" size="icon" aria-label="Choose project folder" title="Choose project folder" onClick={props.onCreateProject}><Icon name="plus"/></Button>
        <Button variant="ghost" size="icon" aria-label="Manage projects" title="Manage projects" onClick={props.onOpenLibrary}><Icon name="folder"/></Button></div>
    </AnimatedSidebarHeader>
    <AnimatedSidebarContent className="workspace-session-list">
      {modules.map(group => <AnimatedSidebarGroup key={group.id}>
        <AnimatedSidebarGroupLabel>{group.label}</AnimatedSidebarGroupLabel>
        {group.children.map(item => <Button key={item.id} variant="ghost" className="product-page-nav"
          disabled={Boolean(item.requiresAnalysis && !props.analysis || item.requiresGsea && !props.gsea)}
          aria-pressed={props.view === item.view && (!item.plotSection || props.plotSection === item.plotSection)}
          onClick={() => { if (item.plotSection) props.onPlotSection(item.plotSection); props.onNavigate(item.view); }}>
          <Icon name={item.icon}/>{item.label}</Button>)}
      </AnimatedSidebarGroup>)}
      {Boolean(props.project.analyses?.length) && <AnimatedSidebarGroup><AnimatedSidebarGroupLabel>{props.project.name}</AnimatedSidebarGroupLabel>
        {props.project.analyses?.map(item => <div key={item.id}>
          <Button variant="ghost" className="product-page-nav" onClick={() => props.onOpenAnalysis(item.id)}><Icon name="flask"/>{item.name}</Button>
          {item.gseaRuns.map(run => <Button key={run.id} variant="ghost" className="product-page-nav nested" onClick={() => props.onOpenGsea(item.id, run.id)}><Icon name="bar"/>{run.name}</Button>)}
        </div>)}
      </AnimatedSidebarGroup>}
    </AnimatedSidebarContent>
    <div className="outer-nav-footer"><Button className="product-page-nav" variant="ghost" aria-pressed={props.view === "settings"} onClick={() => props.onNavigate("settings")}><Icon name="settings"/>Runtime settings</Button></div>
    {props.open && <WorkspaceResizeHandle label="Resize navigation" controls="workspace-project-sidebar" target={() => sidebar.current}
      variable="--sidebar-width" edge="right" className="workspace-resize-sidebar" limits={() => ({min:160,max:260})} onCommit={props.onWidth} onReset={() => props.onWidth(196)}/>}
  </AnimatedSidebar>;
}
