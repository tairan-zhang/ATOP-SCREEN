import {Button} from "@atop/ui";
import { CenterMorphModal, CenterMorphModalContent } from "@atop/ui/motion/center-morph-modal";
import { useEffect, useRef, useState } from "react";
import { desktopApi } from "../services/desktop-api";
import type { ProjectFolderInspection } from "../types/atop";
import { Icon } from "./Icon";

type ProjectDraft = {
  name: string;
  rootPath: string;
};

type ProjectDialogProps = {
  open: boolean;
  onClose: () => void;
  onCreate: (draft: ProjectDraft) => void;
  required?: boolean;
};

export function ProjectDialog({ open, onClose, onCreate, required = false }: ProjectDialogProps) {
  const [name, setName] = useState("");
  const [rootPath, setRootPath] = useState("");
  const [choosingFolder, setChoosingFolder] = useState(false);
  const [inspection, setInspection] = useState<ProjectFolderInspection | null>(null);
  const [error, setError] = useState("");
  const nameInput = useRef<HTMLInputElement>(null);

  useEffect(() => {
    if (!open) return;
    setName("");
    setRootPath("");
    setInspection(null);
    setError("");
  }, [open]);

  async function chooseFolder() {
    setChoosingFolder(true);
    try {
      const selected = await desktopApi.pickDirectory();
      if (selected) {
        setError("");
        const folderInspection = await desktopApi.projects.inspectFolder(selected);
        setRootPath(folderInspection.folderPath);
        setName(folderInspection.name);
        setInspection(folderInspection);
        window.requestAnimationFrame(() => nameInput.current?.select());
      }
    } catch (inspectionError) {
      setInspection(null);
      setError(inspectionError instanceof Error ? inspectionError.message : String(inspectionError));
    } finally {
      setChoosingFolder(false);
    }
  }

  if (!open) return null;

  return (
    <CenterMorphModal open onOpenChange={next => {if (!next && !required) onClose();}}>
      <CenterMorphModalContent ariaLabel="Choose a project folder" className="screen-modal" dismissible={!required} showCloseButton={!required}>
      <form
        aria-labelledby="project-dialog-title"
        className="project-dialog"
        onSubmit={(event) => {
          event.preventDefault();
          const trimmedName = name.trim();
          if (trimmedName && rootPath && inspection) onCreate({ name: trimmedName, rootPath });
        }}
      >
        <header className="project-dialog-heading"><Icon name="folder" size={23}/><h2 id="project-dialog-title">Choose a project folder</h2></header>
        <div className="project-dialog__body">
          <div className="project-folder-field">
            <Button variant="secondary" className="secondary-button" disabled={choosingFolder} onClick={() => void chooseFolder()} type="button">
              <Icon name="folder" />{choosingFolder ? "Checking…" : rootPath ? "Change Folder…" : "Choose Folder…"}
            </Button>
            {rootPath && <small className="project-folder-path">{rootPath}</small>}
          </div>
          {inspection && <label className="field">
            <span className="field-label">Project name</span>
            <input required disabled={inspection.existing} ref={nameInput} value={name} onChange={(event) => setName(event.target.value)} />
          </label>}
          {inspection?.existing && <small>{inspection.analysisCount} saved analyses · {inspection.historyCount} history records</small>}
          {inspection && !inspection.existing && inspection.hasOtherFiles && <small>Project files will be added to this nonempty folder.</small>}
          {error ? <div className="project-folder-inspection is-error"><Icon name="info" size={15} /><div><strong>Folder could not be opened</strong><span>{error}</span></div></div> : null}
        </div>
        <div className="project-dialog__actions">
          {!required ? <Button variant="secondary" className="secondary-button" onClick={onClose} type="button">Cancel</Button> : null}
          <Button className="primary-button" disabled={!name.trim() || !rootPath || !inspection} type="submit"><Icon name={inspection?.existing ? "folder" : "plus"} />{inspection?.existing ? "Open Project" : "Create Project"}</Button>
        </div>
      </form>
      </CenterMorphModalContent>
    </CenterMorphModal>
  );
}

export type { ProjectDraft };
