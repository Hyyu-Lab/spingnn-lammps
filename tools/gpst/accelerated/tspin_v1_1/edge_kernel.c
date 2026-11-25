/* Fused edge arithmetic; reductions remain in NumPy. No reassociation. */
#include <math.h>
#include <stdint.h>

static void poly(double x, int degree, double *p, double *d) {
    p[0] = 1.0; p[1] = x; d[0] = 0.0; d[1] = 1.0;
    for (int k = 1; k < degree; ++k) {
        p[k+1] = ((2*k+1)*x*p[k] - k*p[k-1])/(k+1);
        d[k+1] = ((2*k+1)*(p[k]+x*d[k]) - k*d[k-1])/(k+1);
    }
}

static double dot(const double *a, const double *b) {
    return (a[0]*b[0] + a[1]*b[1]) + a[2]*b[2];
}

static void cross(const double *a, const double *b, double *c) {
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

int tspin_edges(int64_t count, int64_t sites, int64_t edges,
    const double *s, const double *f, const int64_t *left, const int64_t *right,
    const double *w, const double *u, const double *zdir, const double *zforce,
    double *sn, double *sd,
    double *zn, double *zd, double *ph, double *cn, double *ss, double *an) {
    const int orders[4] = {1,2,4,8};
    for (int64_t t = 0; t < count; ++t) {
        for (int64_t e = 0; e < edges; ++e) {
            const double *a = s + (t*sites+left[e])*3;
            const double *b = s + (t*sites+right[e])*3;
            const double *fa = f + (t*sites+left[e])*3;
            const double *fb = f + (t*sites+right[e])*3;
            const double *v = u + e*3;
            double x = dot(a,b), p[9], d[9], cparts[3];
            if (!isfinite(x) || fabs(x) > 1.0+1e-10) return 1;
            x = fmin(1.0, fmax(-1.0,x));
            poly(x,8,p,d);
            for (int k=0; k<3; ++k) cparts[k]=b[k]*fa[k]+a[k]*fb[k];
            double c=(cparts[0]+cparts[1])+cparts[2];
            for (int k=0; k<4; ++k) {
                sn[(t*4+k)*edges+e]=(w[e]*d[orders[k]])*c;
                sd[(t*4+k)*edges+e]=w[e]*p[orders[k]];
            }
            for (int k=0; k<2; ++k) {
                double z=zdir[(t*2+k)*edges+e], pp[9], dd[9];
                if (!isfinite(z) || fabs(z) > 1.0+1e-10) return 1;
                z=fmin(1.0,fmax(-1.0,z));
                poly(z,4,pp,dd);
                zn[(k*edges+e)*count+t]=(w[e]*dd[4])*zforce[(t*2+k)*edges+e];
                zd[(k*edges+e)*count+t]=w[e]*pp[4];
            }
            double ab[3], bv[3], va[3];
            cross(a,b,ab); cross(b,v,bv); cross(v,a,va);
            ph[t*edges+e]=w[e]*dot(v,ab);
            cn[t*edges+e]=w[e]*(dot(fa,bv)+dot(fb,va));
            for (int k=0; k<3; ++k) for (int l=0; l<3; ++l) {
                int64_t q=((t*edges+e)*3+k)*3+l;
                ss[q]=(w[e]*a[k])*b[l];
                an[q]=w[e]*(fa[k]*b[l]+a[k]*fb[l]);
            }
        }
    }
    return 0;
}
