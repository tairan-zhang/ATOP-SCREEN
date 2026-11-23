import {Input} from "@atop/ui";
import {Checkbox} from "@atop/ui/motion/checkbox";

export function TextField({ label, hint, value, onChange }: { label: string; hint?: string; value: string; onChange: (value: string) => void }) {
  return <div className="field"><Input label={label} value={value} onChange={onChange}/>{hint && <small>{hint}</small>}</div>;
}

export function CheckboxField({ label, hint, checked, disabled = false, onChange }: { label: string; hint?: string; checked: boolean; disabled?: boolean; onChange: (value: boolean) => void }) {
  return <div className="field"><Checkbox label={label} checked={checked} disabled={disabled} onCheckedChange={onChange}/>{hint && <small>{hint}</small>}</div>;
}

export function NumberField({ label, hint, value, step = "1", min = 0, max, onChange }: { label: string; hint?: string; value: number; step?: string; min?: number; max?: number; onChange: (value: number) => void }) {
  return <div className="field"><Input label={label} aria-label={label} type="number" required min={min} max={max} step={step}
    value={Number.isFinite(value) ? String(value) : ""} onChange={text => onChange(text === "" ? NaN : Number(text))}/>
    {hint && <small>{hint}</small>}</div>;
}

export { SelectField } from "./SelectField";
export { MultiSelectField } from "./MultiSelectField";
export { ColorField } from "./ColorField";
