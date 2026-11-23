# ATOP SCREEN icon

The macOS application uses `ATOP_black.icns`, the black, white and gray Liquid Glass artwork from `ATOP_black.icon`. The Icon Composer document contains the editable glass layers.

The sidebar and page icon use the 128-pixel rendition from the same icon set. To refresh it from the software root on macOS:

```bash
icon_build=$(mktemp -d)
iconutil --convert iconset --output "$icon_build/ATOP.iconset" packaging/icons/ATOP_black.icns
cp "$icon_build/ATOP.iconset/icon_128x128.png" web/public/brand/atop-glass.png
```
