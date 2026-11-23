import { Select, SelectTrigger, SelectContent, SelectItem } from "@atop/ui/motion/select";

export type SelectFieldProps = {
  label: string;
  value: string;
  options: string[];
  optionLabels?: Record<string, string>;
  placeholder?: string;
  hint?: string;
  requirement?: "Required" | "Optional";
  disabled?: boolean;
  onChange: (value: string) => void;
};

export function SelectField({label, options, optionLabels, placeholder, hint, value, onChange, disabled}: SelectFieldProps) {
  return <div className="field" role="group" aria-label={label}><span>{label}</span>
    <Select value={value} onValueChange={onChange} disabled={disabled || !options.length}>
      <SelectTrigger>{optionLabels?.[value] ?? (value || placeholder || "Choose…")}</SelectTrigger>
      <SelectContent>{options.map(option => <SelectItem key={option} value={option}>{optionLabels?.[option] ?? option}</SelectItem>)}</SelectContent>
    </Select>{hint && <small>{hint}</small>}
  </div>;
}
