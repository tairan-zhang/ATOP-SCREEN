import { useEffect, useState } from "react";

type ColorFieldProps = {
  hint?: string;
  label: string;
  onChange: (value: string) => void;
  value: string;
};

function normalizeHexColor(value: string) {
  const compact = value.trim();
  const short = compact.match(/^#?([0-9a-f]{3})$/i);
  if (short) return `#${short[1].split("").map((part) => `${part}${part}`).join("").toUpperCase()}`;
  const full = compact.match(/^#?([0-9a-f]{6})$/i);
  return full ? `#${full[1].toUpperCase()}` : null;
}

export function ColorField({ hint, label, onChange, value }: ColorFieldProps) {
  const [draft, setDraft] = useState(value.toUpperCase());
  useEffect(() => setDraft(value.toUpperCase()), [value]);

  const commitDraft = () => {
    const normalized = normalizeHexColor(draft);
    if (normalized) {
      onChange(normalized);
      setDraft(normalized);
    } else {
      setDraft(value.toUpperCase());
    }
  };

  return (
    <div className="field color-field" role="group" aria-label={label}>
      <span className="field-label">{label}</span>
      <span className="color-control">
        <input aria-label={`${label} picker`} onChange={(event) => onChange(event.target.value.toUpperCase())} type="color" value={value} />
        <input
          aria-label={`${label} value`}
          className="color-control__value"
          maxLength={7}
          onBlur={commitDraft}
          onChange={(event) => setDraft(event.target.value)}
          onKeyDown={(event) => { if (event.key === "Enter") event.currentTarget.blur(); }}
          spellCheck={false}
          value={draft}
        />
      </span>
      {hint && <small>{hint}</small>}
    </div>
  );
}
