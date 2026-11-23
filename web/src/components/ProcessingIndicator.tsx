export function ProcessingIndicator({label}: {label:string}) {
  return <span role="status" aria-live="polite">{label}</span>;
}
