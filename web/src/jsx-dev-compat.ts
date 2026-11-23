// atop-ui 0.1.0 ships JSX development calls in its distributed modules.
// Resolve that runtime to React's production JSX without changing UI components.
import { jsx, jsxs, Fragment } from 'react/jsx-runtime';
export { jsx, jsxs, Fragment };
export function jsxDEV(type: Parameters<typeof jsx>[0], props: Parameters<typeof jsx>[1], key?: string, staticChildren?: boolean) {
 return (staticChildren ? jsxs : jsx)(type, props, key);
}
