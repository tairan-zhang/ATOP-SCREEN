import {Button} from "@atop/ui";
import { useState } from "react";
import type { DragEvent } from "react";
import { Icon } from "./Icon";
import { buildReplicatePairs, reorderReplicates, validateReplicatePairing } from "../features/analysis/replicatePairing";

type Group = "treatment" | "control";

type Props = {
  control: string[];
  disabled?: boolean;
  onControlChange: (values: string[]) => void;
  onTreatmentChange: (values: string[]) => void;
  treatment: string[];
};

const pairStatusCopy = {
  matched: { icon: "check", label: "Replicate number matches" },
  review: { icon: "info", label: "Review manually" },
  mismatch: { icon: "info", label: "Check replicate order" },
  missing: { icon: "info", label: "Missing pair" },
} as const;

export function ReplicatePairingEditor({ control, disabled = false, onControlChange, onTreatmentChange, treatment }: Props) {
  const [dragged, setDragged] = useState<{ group: Group; index: number } | null>(null);
  const [dropTarget, setDropTarget] = useState<{ group: Group; index: number } | null>(null);
  const validation = validateReplicatePairing(treatment, control);
  const pairs = buildReplicatePairs(treatment, control);

  const valuesFor = (group: Group) => group === "treatment" ? treatment : control;
  const updateGroup = (group: Group, values: string[]) => group === "treatment" ? onTreatmentChange(values) : onControlChange(values);

  const move = (group: Group, fromIndex: number, toIndex: number) => {
    updateGroup(group, reorderReplicates(valuesFor(group), fromIndex, toIndex));
  };

  const drop = (event: DragEvent<HTMLElement>, group: Group, toIndex: number) => {
    event.preventDefault();
    if (dragged?.group === group) move(group, dragged.index, toIndex);
    setDragged(null);
    setDropTarget(null);
  };

  function renderSample(group: Group, index: number, column: string | null) {
    if (!column) return <span className="pair-missing">Missing</span>;
    const values = valuesFor(group);
    return <div className="pair-sample" draggable={!disabled}
      data-drop-target={dropTarget?.group === group && dropTarget.index === index}
      onDragStart={event => {
        event.dataTransfer.effectAllowed = "move";
        event.dataTransfer.setData("text/plain", column);
        setDragged({group, index});
      }}
      onDragOver={event => {
        if (dragged?.group === group) { event.preventDefault(); setDropTarget({group, index}); }
      }}
      onDrop={event => drop(event, group, index)}
      onDragEnd={() => { setDragged(null); setDropTarget(null); }}>
      <span className="pair-name" title={column}>{column}</span>
      <div className="pair-actions">
        <Button variant="ghost" type="button" aria-label={`Move ${column} up`} disabled={disabled || index === 0} onClick={() => move(group,index,index-1)}><Icon name="arrowUp" size={12}/></Button>
        <Button variant="ghost" type="button" aria-label={`Move ${column} down`} disabled={disabled || index === values.length-1} onClick={() => move(group,index,index+1)}><Icon name="chevronDown" size={12}/></Button>
        <Button variant="ghost" type="button" aria-label={`Remove ${column}`} disabled={disabled} onClick={() => updateGroup(group,values.filter(v => v !== column))}><Icon name="close" size={12}/></Button>
      </div>
    </div>;
  }

  if (!pairs.length) return null;
  return <div className="pair-editor">
    <div className="pair-grid pair-heading"><span>Endpoint order</span><span>Baseline order</span><span/></div>
    {pairs.map(pair => <div className="pair-grid" key={pair.index}>
      {renderSample("treatment",pair.index,pair.treatment)}
      {renderSample("control",pair.index,pair.control)}
      <span className={`pair-check is-${pair.status}`} title={pairStatusCopy[pair.status].label} aria-label={pairStatusCopy[pair.status].label}>
        <Icon name={pairStatusCopy[pair.status].icon} size={14}/>
      </span>
    </div>)}
    {!validation.ready && <small role="status">{validation.overlap.length ? "A column is used in both groups." : "Matched pairs need equal numbers of endpoints and baselines."}</small>}
  </div>;
}
