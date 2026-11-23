import { Children, isValidElement, type ReactElement, type ReactNode } from "react";
import { describe, expect, it, vi } from "vitest";
import { WorkspaceTitlebar } from "./WorkspaceTitlebar";

type ButtonProps = { "aria-label"?: string; className?: string; onClick?: () => void; children?: ReactNode };

function buttonsIn(node: ReactNode): ReactElement<ButtonProps>[] {
  const buttons: ReactElement<ButtonProps>[] = [];
  Children.forEach(node, (child) => {
    if (!isValidElement<ButtonProps>(child)) return;
    if (child.props["aria-label"]) buttons.push(child);
    buttons.push(...buttonsIn(child.props.children));
  });
  return buttons;
}

describe("WorkspaceTitlebar", () => {
  it("keeps both pane toggles clickable when the sidebar and preview are collapsed", () => {
    const onToggleSidebar = vi.fn();
    const onTogglePreview = vi.fn();
    const titlebar = WorkspaceTitlebar({
      onTogglePreview,
      onToggleSidebar,
      pageTitle: "Screen analysis",
      previewOpen: false,
      projectName: "Test project",
      sidebarOpen: false,
    });
    const buttons = buttonsIn(titlebar);

    expect(buttons.map((button) => button.props["aria-label"])).toEqual(["Show navigation", "Show preview"]);
    buttons[0].props.onClick?.();
    buttons[1].props.onClick?.();
    expect(onToggleSidebar).toHaveBeenCalledOnce();
    expect(onTogglePreview).toHaveBeenCalledOnce();
  });
});
