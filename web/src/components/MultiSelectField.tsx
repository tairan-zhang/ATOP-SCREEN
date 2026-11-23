import {MultiSelect, MultiSelectContent, MultiSelectTrigger, MultiSelectValue, MultiSelectItem, MultiSelectList} from "@atop/ui/motion/multi-select";

type MultiSelectFieldProps = {
  label: string;
  values: string[];
  options: string[];
  placeholder?: string;
  hint?: string;
  requirement?: "Required" | "Optional";
  disabled?: boolean;
  onChange: (values: string[]) => void;
};

export function MultiSelectField({label, values, options, placeholder, hint, disabled, onChange}: MultiSelectFieldProps) {
  return <div className="field" role="group" aria-label={label}><span>{label}</span>
    <MultiSelect value={values} onValueChange={onChange} disabled={disabled || !options.length}>
      <MultiSelectTrigger><MultiSelectValue placeholder={placeholder ?? "Choose columns…"}/></MultiSelectTrigger>
      <MultiSelectContent><MultiSelectList ariaLabel={label}>{options.map(value => <MultiSelectItem key={value} value={value}>{value}</MultiSelectItem>)}</MultiSelectList></MultiSelectContent>
    </MultiSelect>{hint && <small>{hint}</small>}
  </div>;
}
