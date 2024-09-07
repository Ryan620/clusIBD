import { EllipseOval, EllipseOvalView } from "./ellipse_oval";
export class OvalView extends EllipseOvalView {
    _map_data() {
        let sw;
        const n = this._x.length;
        this.sw = new Float64Array(n);
        if (this.model.properties.width.units == "data")
            sw = this.sdist(this.renderer.xscale, this._x, this._width, 'center');
        else
            sw = this._width;
        // oval drawn from bezier curves = ellipse with width reduced by 3/4
        for (let i = 0; i < n; i++)
            this.sw[i] = sw[i] * 0.75;
        if (this.model.properties.height.units == "data")
            this.sh = this.sdist(this.renderer.yscale, this._y, this._height, 'center');
        else
            this.sh = this._height;
    }
}
OvalView.__name__ = "OvalView";
export class Oval extends EllipseOval {
    constructor(attrs) {
        super(attrs);
    }
    static init_Oval() {
        this.prototype.default_view = OvalView;
    }
}
Oval.__name__ = "Oval";
Oval.init_Oval();
//# sourceMappingURL=oval.js.map