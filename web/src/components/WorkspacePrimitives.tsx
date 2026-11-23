import { Icon } from "./Icon";

export function Prerequisite({ children }: { children: string }) {
  return <div className="prerequisite"><Icon name="info" size={16} /><span>{children}</span></div>;
}

export function ErrorBanner({ message }: { message: string }) {
  return <div className="error-banner"><Icon name="info" /><span>{message}</span></div>;
}
