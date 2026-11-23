import type {ReactElement, ReactNode} from "react";
import {Tooltip} from "@atop/ui/motion/tooltip";
export function Tip({label, children, disabled}: {label:ReactNode; children:ReactElement; disabled?:boolean}) {
  return disabled || !label ? <>{children}</> : <Tooltip content={label} side="top">{children}</Tooltip>;
}
