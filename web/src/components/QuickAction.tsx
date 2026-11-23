import {Button} from "@atop/ui";
import { Icon } from "./Icon";
import { Tip } from "./Tooltip";

type QuickActionProps = {
  icon: string;
  label: string;
  onClick: () => void;
  disabled?: boolean;
  disabledReason?: string;
};

export function QuickAction({ icon, label, onClick, disabled = false, disabledReason = "" }: QuickActionProps) {
  return (
    <Tip label={disabledReason} disabled={!disabled || !disabledReason}>
      <span className="quick-action-tip">
        <Button variant="secondary" className="quick-action" disabled={disabled} aria-disabled={disabled} onClick={onClick} type="button">
          <span className="quick-action-icon"><Icon name={icon} /></span>
          <span>{label}</span>
          <Icon name="chevronRight" size={15} />
        </Button>
      </span>
    </Tip>
  );
}
