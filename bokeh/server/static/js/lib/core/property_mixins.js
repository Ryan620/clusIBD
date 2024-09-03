import * as p from "./properties";
import { extend } from "./util/object";
function _gen_mixin(mixin, prefix) {
    const result = {};
    for (const name in mixin) {
        const prop = mixin[name];
        result[prefix + name] = prop;
    }
    return result;
}
const _line_mixin = {
    line_color: [p.ColorSpec, 'black'],
    line_width: [p.NumberSpec, 1],
    line_alpha: [p.NumberSpec, 1.0],
    line_join: [p.LineJoin, 'bevel'],
    line_cap: [p.LineCap, 'butt'],
    line_dash: [p.Array, []],
    line_dash_offset: [p.Number, 0],
};
export const line = (prefix = "") => _gen_mixin(_line_mixin, prefix);
const _fill_mixin = {
    fill_color: [p.ColorSpec, 'gray'],
    fill_alpha: [p.NumberSpec, 1.0],
};
export const fill = (prefix = "") => _gen_mixin(_fill_mixin, prefix);
const _hatch_mixin = {
    hatch_color: [p.ColorSpec, 'black'],
    hatch_alpha: [p.NumberSpec, 1.0],
    hatch_scale: [p.NumberSpec, 12.0],
    hatch_pattern: [p.StringSpec, null],
    hatch_weight: [p.NumberSpec, 1.0],
    hatch_extra: [p.Any, {}],
};
export const hatch = (prefix = "") => _gen_mixin(_hatch_mixin, prefix);
const _text_mixin = {
    text_font: [p.Font, 'helvetica'],
    text_font_size: [p.FontSizeSpec, '12pt'],
    text_font_style: [p.FontStyle, 'normal'],
    text_color: [p.ColorSpec, '#444444'],
    text_alpha: [p.NumberSpec, 1.0],
    text_align: [p.TextAlign, 'left'],
    text_baseline: [p.TextBaseline, 'bottom'],
    text_line_height: [p.Number, 1.2],
};
export const text = (prefix = "") => _gen_mixin(_text_mixin, prefix);
export function create(configs) {
    const result = {};
    for (const config of configs) {
        const [kind, prefix] = config.split(":");
        let mixin;
        switch (kind) {
            case "line":
                mixin = line;
                break;
            case "fill":
                mixin = fill;
                break;
            case "hatch":
                mixin = hatch;
                break;
            case "text":
                mixin = text;
                break;
            default:
                throw new Error(`Unknown property mixin kind '${kind}'`);
        }
        extend(result, mixin(prefix));
    }
    return result;
}
//# sourceMappingURL=property_mixins.js.map