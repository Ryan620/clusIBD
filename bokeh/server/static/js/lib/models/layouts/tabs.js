import { Grid, ContentBox } from "../../core/layout";
import { div, position, size, scroll_size, show, hide, display, undisplay, children } from "../../core/dom";
import { sum, remove_at } from "../../core/util/array";
import * as p from "../../core/properties";
import { LayoutDOM, LayoutDOMView } from "./layout_dom";
import { Model } from "../../model";
import { bk_left, bk_right, bk_active, bk_side } from "../../styles/mixins";
import { bk_tabs_header, bk_headers, bk_headers_wrapper, bk_tab, bk_close } from "../../styles/tabs";
import { bk_btn, bk_btn_default, bk_btn_group } from "../../styles/buttons";
import { bk_caret } from "../../styles/menus";
export class TabsView extends LayoutDOMView {
    connect_signals() {
        super.connect_signals();
        this.connect(this.model.properties.tabs.change, () => this.rebuild());
        this.connect(this.model.properties.active.change, () => this.on_active_change());
    }
    get child_models() {
        return this.model.tabs.map((tab) => tab.child);
    }
    _update_layout() {
        const loc = this.model.tabs_location;
        const vertical = loc == "above" || loc == "below";
        // XXX: this is a hack, this should be handled by "fit" policy in grid layout
        const { scroll_el, headers_el } = this;
        this.header = new class extends ContentBox {
            _measure(viewport) {
                const min_headers = 3;
                const scroll = size(scroll_el);
                const headers = children(headers_el).slice(0, min_headers).map((el) => size(el));
                const { width, height } = super._measure(viewport);
                if (vertical) {
                    const min_width = scroll.width + sum(headers.map((size) => size.width));
                    return { width: viewport.width != Infinity ? viewport.width : min_width, height };
                }
                else {
                    const min_height = scroll.height + sum(headers.map((size) => size.height));
                    return { width, height: viewport.height != Infinity ? viewport.height : min_height };
                }
            }
        }(this.header_el);
        if (vertical)
            this.header.set_sizing({ width_policy: "fit", height_policy: "fixed" });
        else
            this.header.set_sizing({ width_policy: "fixed", height_policy: "fit" });
        let row = 1;
        let col = 1;
        switch (loc) {
            case "above":
                row -= 1;
                break;
            case "below":
                row += 1;
                break;
            case "left":
                col -= 1;
                break;
            case "right":
                col += 1;
                break;
        }
        const header = { layout: this.header, row, col };
        const panels = this.child_views.map((child_view) => {
            return { layout: child_view.layout, row: 1, col: 1 };
        });
        this.layout = new Grid([header, ...panels]);
        this.layout.set_sizing(this.box_sizing());
    }
    update_position() {
        super.update_position();
        this.header_el.style.position = "absolute"; // XXX: do it in position()
        position(this.header_el, this.header.bbox);
        const loc = this.model.tabs_location;
        const vertical = loc == "above" || loc == "below";
        const scroll_el_size = size(this.scroll_el);
        const headers_el_size = scroll_size(this.headers_el);
        if (vertical) {
            const { width } = this.header.bbox;
            if (headers_el_size.width > width) {
                this.wrapper_el.style.maxWidth = `${width - scroll_el_size.width}px`;
                display(this.scroll_el);
            }
            else {
                this.wrapper_el.style.maxWidth = "";
                undisplay(this.scroll_el);
            }
        }
        else {
            const { height } = this.header.bbox;
            if (headers_el_size.height > height) {
                this.wrapper_el.style.maxHeight = `${height - scroll_el_size.height}px`;
                display(this.scroll_el);
            }
            else {
                this.wrapper_el.style.maxHeight = "";
                undisplay(this.scroll_el);
            }
        }
        const { child_views } = this;
        for (const child_view of child_views)
            hide(child_view.el);
        const tab = child_views[this.model.active];
        if (tab != null)
            show(tab.el);
    }
    render() {
        super.render();
        const { active } = this.model;
        const loc = this.model.tabs_location;
        const vertical = loc == "above" || loc == "below";
        const headers = this.model.tabs.map((tab, i) => {
            const el = div({ class: [bk_tab, i == active ? bk_active : null] }, tab.title);
            el.addEventListener("click", (event) => {
                if (event.target == event.currentTarget)
                    this.change_active(i);
            });
            if (tab.closable) {
                const close_el = div({ class: bk_close });
                close_el.addEventListener("click", (event) => {
                    if (event.target == event.currentTarget) {
                        this.model.tabs = remove_at(this.model.tabs, i);
                        const ntabs = this.model.tabs.length;
                        if (this.model.active > ntabs - 1)
                            this.model.active = ntabs - 1;
                    }
                });
                el.appendChild(close_el);
            }
            return el;
        });
        this.headers_el = div({ class: [bk_headers] }, headers);
        this.wrapper_el = div({ class: bk_headers_wrapper }, this.headers_el);
        const left_el = div({ class: [bk_btn, bk_btn_default], disabled: "" }, div({ class: [bk_caret, bk_left] }));
        const right_el = div({ class: [bk_btn, bk_btn_default] }, div({ class: [bk_caret, bk_right] }));
        let scroll_index = 0;
        const do_scroll = (dir) => {
            return () => {
                const ntabs = this.model.tabs.length;
                if (dir == "left")
                    scroll_index = Math.max(scroll_index - 1, 0);
                else
                    scroll_index = Math.min(scroll_index + 1, ntabs - 1);
                if (scroll_index == 0)
                    left_el.setAttribute("disabled", "");
                else
                    left_el.removeAttribute("disabled");
                if (scroll_index == ntabs - 1)
                    right_el.setAttribute("disabled", "");
                else
                    right_el.removeAttribute("disabled");
                const sizes = children(this.headers_el)
                    .slice(0, scroll_index)
                    .map((el) => el.getBoundingClientRect());
                if (vertical) {
                    const left = -sum(sizes.map((size) => size.width));
                    this.headers_el.style.left = `${left}px`;
                }
                else {
                    const top = -sum(sizes.map((size) => size.height));
                    this.headers_el.style.top = `${top}px`;
                }
            };
        };
        left_el.addEventListener("click", do_scroll("left"));
        right_el.addEventListener("click", do_scroll("right"));
        this.scroll_el = div({ class: bk_btn_group }, left_el, right_el);
        this.header_el = div({ class: [bk_tabs_header, bk_side(loc)] }, this.scroll_el, this.wrapper_el);
        this.el.appendChild(this.header_el);
    }
    change_active(i) {
        if (i != this.model.active) {
            this.model.active = i;
        }
    }
    on_active_change() {
        const i = this.model.active;
        const headers = children(this.headers_el);
        for (const el of headers)
            el.classList.remove(bk_active);
        headers[i].classList.add(bk_active);
        const { child_views } = this;
        for (const child_view of child_views)
            hide(child_view.el);
        show(child_views[i].el);
    }
}
TabsView.__name__ = "TabsView";
export class Tabs extends LayoutDOM {
    constructor(attrs) {
        super(attrs);
    }
    static init_Tabs() {
        this.prototype.default_view = TabsView;
        this.define({
            tabs: [p.Array, []],
            tabs_location: [p.Location, "above"],
            active: [p.Number, 0],
        });
    }
}
Tabs.__name__ = "Tabs";
Tabs.init_Tabs();
export class Panel extends Model {
    constructor(attrs) {
        super(attrs);
    }
    static init_Panel() {
        this.define({
            title: [p.String, ""],
            child: [p.Instance],
            closable: [p.Boolean, false],
        });
    }
}
Panel.__name__ = "Panel";
Panel.init_Panel();
//# sourceMappingURL=tabs.js.map