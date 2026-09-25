// Independent probe harness for the isolated Taylor-jet implementation.
#define main direct_study_main
#include "direct_derivatives.cpp"
#undef main
int main() {
    std::cout << std::setprecision(17);
    for (bool coupled : {false, true}) {
        direct::Geometry g{.12,.18,.20,.06,-.03,.22,coupled?1.2:.7,coupled};
        double residual=0, h=.25;
        auto rs=g.InverseJet(.23,.31,h,residual);
        // A tensor polynomial with nontrivial mixed derivatives through order 3.
        auto r=rs[0],s=rs[1];auto u=r*r*r*s + .4*r*s*s + .7*r*r - .2*s;
        for(int i=1;i<direct::count;++i)
            std::cout<<(coupled?"coupled":"separable")<<','<<direct::ax[i]<<','<<direct::ay[i]<<','
                     <<u.c[i]*direct::Factorial(direct::ax[i])*direct::Factorial(direct::ay[i])/std::pow(h,direct::ax[i]+direct::ay[i])<<'\n';
    }
}
