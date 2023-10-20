#
# This is an autogenerated file! It was created on 2023-10-19.
# 

const expressions = Base.ImmutableDict(
	"HarmonicOscillatorPotential" => L"\left( 0.5 \omega^{2} x^{2}, \  \left\{ \mathtt{\text{x}} : x\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{omega}} : \omega\right\}\right)",
	"HenonHeilesPotential" => L"\left( 1.0 x^{2} y + 0.5 x^{2} - 0.333333333333333 y^{3} + 0.5 y^{2}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y\right\}, \  \left\{ \mathtt{\text{G}} : G\right\}\right)",
	"HernquistPotential" => L"\left( - \frac{G m}{c + \sqrt{x^{2} + y^{2} + z^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{c}} : c, \  \mathtt{\text{m}} : m\right\}\right)",
	"IsochronePotential" => L"\left( - \frac{G m}{b + \sqrt{b^{2} + x^{2} + y^{2} + z^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{b}} : b, \  \mathtt{\text{m}} : m\right\}\right)",
	"JaffePotential" => L"\left( \frac{G m \log{\left(\frac{\sqrt{x^{2} + y^{2} + z^{2}}}{c + \sqrt{x^{2} + y^{2} + z^{2}}} \right)}}{c}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{c}} : c, \  \mathtt{\text{m}} : m\right\}\right)",
	"KeplerPotential" => L"\left( - \frac{G m}{\sqrt{x^{2} + y^{2} + z^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{m}} : m\right\}\right)",
	"KuzminPotential" => L"\left( - \frac{G m}{\sqrt{x^{2} + y^{2} + \left(a + \left|{z}\right|\right)^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{a}} : a, \  \mathtt{\text{m}} : m\right\}\right)",
	"LogarithmicPotential" => L"\left( 0.5 v_{c}^{2} \log{\left(r_{h}^{2} + \frac{z^{2}}{q_{3}^{2}} + \frac{y^{2}}{q_{2}^{2}} + \frac{x^{2}}{q_{1}^{2}} \right)}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{phi}} : \phi, \  \mathtt{\text{q1}} : q_{1}, \  \mathtt{\text{q2}} : q_{2}, \  \mathtt{\text{q3}} : q_{3}, \  \mathtt{\text{r\_h}} : r_{h}, \  \mathtt{\text{v\_c}} : v_{c}\right\}\right)",
	"LongMuraliBarPotential" => L"\left( \frac{G m \log{\left(\frac{- a + x \cos{\left(\alpha \right)} + y \sin{\left(\alpha \right)} + \sqrt{\left(b + \sqrt{c^{2} + z^{2}}\right)^{2} + \left(- x \sin{\left(\alpha \right)} + y \cos{\left(\alpha \right)}\right)^{2} + \left(a - x \cos{\left(\alpha \right)} - y \sin{\left(\alpha \right)}\right)^{2}}}{a + x \cos{\left(\alpha \right)} + y \sin{\left(\alpha \right)} + \sqrt{\left(b + \sqrt{c^{2} + z^{2}}\right)^{2} + \left(- x \sin{\left(\alpha \right)} + y \cos{\left(\alpha \right)}\right)^{2} + \left(a + x \cos{\left(\alpha \right)} + y \sin{\left(\alpha \right)}\right)^{2}}} \right)}}{2 a}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{a}} : a, \  \mathtt{\text{alpha}} : \alpha, \  \mathtt{\text{b}} : b, \  \mathtt{\text{c}} : c, \  \mathtt{\text{m}} : m\right\}\right)",
	"MiyamotoNagaiPotential" => L"\left( - \frac{G m}{\sqrt{x^{2} + y^{2} + \left(a + \sqrt{b^{2} + z^{2}}\right)^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{a}} : a, \  \mathtt{\text{b}} : b, \  \mathtt{\text{m}} : m\right\}\right)",
	"NFWPotential" => L"\left( - \frac{G m \log{\left(1 + \frac{\sqrt{\frac{z^{2}}{c^{2}} + \frac{y^{2}}{b^{2}} + \frac{x^{2}}{a^{2}}}}{r_{s}} \right)}}{\sqrt{\frac{z^{2}}{c^{2}} + \frac{y^{2}}{b^{2}} + \frac{x^{2}}{a^{2}}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{a}} : a, \  \mathtt{\text{b}} : b, \  \mathtt{\text{c}} : c, \  \mathtt{\text{m}} : m, \  \mathtt{\text{r\_s}} : r_{s}\right\}\right)",
	"PlummerPotential" => L"\left( - \frac{G m}{\sqrt{b^{2} + x^{2} + y^{2} + z^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{b}} : b, \  \mathtt{\text{m}} : m\right\}\right)",
	"PowerLawCutoffPotential" => L"\left( \frac{G \alpha m \gamma\left(1.5 - \frac{\alpha}{2}, \frac{x^{2} + y^{2} + z^{2}}{r_{c}^{2}}\right)}{2 \sqrt{x^{2} + y^{2} + z^{2}} \Gamma\left(2.5 - \frac{\alpha}{2}\right)} - \frac{3 G m \gamma\left(1.5 - \frac{\alpha}{2}, \frac{x^{2} + y^{2} + z^{2}}{r_{c}^{2}}\right)}{2 \sqrt{x^{2} + y^{2} + z^{2}} \Gamma\left(2.5 - \frac{\alpha}{2}\right)} + \frac{G m \gamma\left(1 - \frac{\alpha}{2}, \frac{x^{2} + y^{2} + z^{2}}{r_{c}^{2}}\right)}{r_{c} \Gamma\left(1.5 - \frac{\alpha}{2}\right)}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{alpha}} : \alpha, \  \mathtt{\text{m}} : m, \  \mathtt{\text{r\_c}} : r_{c}\right\}\right)",
	"SatohPotential" => L"\left( - \frac{G m}{\sqrt{a \left(a + 2 \sqrt{b^{2} + z^{2}}\right) + x^{2} + y^{2} + z^{2}}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{a}} : a, \  \mathtt{\text{b}} : b, \  \mathtt{\text{m}} : m\right\}\right)",
	"StonePotential" => L"\left( - \frac{2 G m \left(- \frac{r_{c} \operatorname{atan}{\left(\frac{\sqrt{x^{2} + y^{2} + z^{2}}}{r_{c}} \right)}}{\sqrt{x^{2} + y^{2} + z^{2}}} + \frac{r_{h} \operatorname{atan}{\left(\frac{\sqrt{x^{2} + y^{2} + z^{2}}}{r_{h}} \right)}}{\sqrt{x^{2} + y^{2} + z^{2}}} + 0.5 \log{\left(\frac{r_{h}^{2} + x^{2} + y^{2} + z^{2}}{r_{c}^{2} + x^{2} + y^{2} + z^{2}} \right)}\right)}{- 3.14159265358979 r_{c} + 3.14159265358979 r_{h}}, \  \left\{ \mathtt{\text{x}} : x, \  \mathtt{\text{y}} : y, \  \mathtt{\text{z}} : z\right\}, \  \left\{ \mathtt{\text{G}} : G, \  \mathtt{\text{m}} : m, \  \mathtt{\text{r\_c}} : r_{c}, \  \mathtt{\text{r\_h}} : r_{h}\right\}\right)",
)
