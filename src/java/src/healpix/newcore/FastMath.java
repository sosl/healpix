package healpix.newcore;

/* Subset of the SLEEF Java library by Naoki Shibata. The code in this file is
   released to the Public Domain. */

/**
 * FastMath class is a Java implementation of the
 * <a href="http://freecode.com/projects/sleef">SLEEF</a>
 * library. Some of the methods can be used as substitutions of the
 * corresponding methods in Math class. They have slightly less
 * accuracy, and some methods are faster compared to those methods in
 * Math class. Please note that the methods in the standard Math class
 * are JNI methods, and the SLEEF library is specialized for SIMD
 * operations.
 */
public class FastMath {
    static double upper(double d) {
        long l = Double.doubleToRawLongBits(d);
        return Double.longBitsToDouble(l & 0xfffffffff8000000L);
    }

    static double mla(double x, double y, double z) { return x * y + z; }

    static double mulsign(double x, double y) { return Math.copySign(1, y) * x; }

    //

    /**
       Returns the absolute value of the argument
     */
    public static double fabs(double d) { return Math.copySign(d, 1); }

    /**
       Checks if the argument is a NaN or not.
     */
    public static boolean isnan(double d) { return d != d; }

    /**
       Checks if the argument is either positive infinity or negative infinity.
     */
    public static boolean isinf(double d) { return fabs(d) == Double.POSITIVE_INFINITY; }

    static boolean ispinf(double d) { return d == Double.POSITIVE_INFINITY; }
    static boolean isminf(double d) { return d == Double.NEGATIVE_INFINITY; }

    /**
       Returns the integer value that is closest to the argument. The
       result is undefined if a denormal number is given.
     */
    public static int rint(double x) { return x < 0 ? (int)(x - 0.5) : (int)(x + 0.5); }

    static double sign(double d) { return Math.copySign(1, d); }

    static double atan2k(double y, double x) {
        double s, t, u;
        int q = 0;

        if (x < 0) { x = -x; q = -2; }
        if (y > x) { t = x; x = y; y = -t; q += 1; }

        s = y / x;
        t = s * s;

        u = -1.88796008463073496563746e-05;
        u = u * t + (0.000209850076645816976906797);
        u = u * t + (-0.00110611831486672482563471);
        u = u * t + (0.00370026744188713119232403);
        u = u * t + (-0.00889896195887655491740809);
        u = u * t + (0.016599329773529201970117);
        u = u * t + (-0.0254517624932312641616861);
        u = u * t + (0.0337852580001353069993897);
        u = u * t + (-0.0407629191276836500001934);
        u = u * t + (0.0466667150077840625632675);
        u = u * t + (-0.0523674852303482457616113);
        u = u * t + (0.0587666392926673580854313);
        u = u * t + (-0.0666573579361080525984562);
        u = u * t + (0.0769219538311769618355029);
        u = u * t + (-0.090908995008245008229153);
        u = u * t + (0.111111105648261418443745);
        u = u * t + (-0.14285714266771329383765);
        u = u * t + (0.199999999996591265594148);
        u = u * t + (-0.333333333333311110369124);

        t = u * t * s + s;
        t = q * (Math.PI/2) + t;

        return t;
    }

    /**
       This method calculates the arc tangent of y/x in radians, using
       the signs of the two arguments to determine the quadrant of the
       result. The results may have maximum error of 2 ulps.
     */
    public static double atan2(double y, double x) {
        double r = atan2k(fabs(y), x);

        r = mulsign(r, x);
        if (isinf(x) || x == 0) r = Math.PI/2 - (isinf(x) ? (sign(x) * (Math.PI  /2)) : 0);
        if (isinf(y)          ) r = Math.PI/2 - (isinf(x) ? (sign(x) * (Math.PI*1/4)) : 0);
        if (            y == 0) r = (sign(x) == -1 ? Math.PI : 0);

        return isnan(x) || isnan(y) ? Double.NaN : mulsign(r, y);
    }

    /**
       This method calculates the arc sine of x in radians. The return
       value is in the range [-pi/2, pi/2]. The results may have
       maximum error of 3 ulps.
     */
    public static double asin(double d) {
        return mulsign(atan2k(fabs(d), Math.sqrt((1+d)*(1-d))), d);
    }

    /**
       This method calculates the arc cosine of x in radians. The
       return value is in the range [0, pi]. The results may have
       maximum error of 3 ulps.
    */
    public static double acos(double d) {
        return mulsign(atan2k(Math.sqrt((1+d)*(1-d)), fabs(d)), d) + (d < 0 ? Math.PI : 0);
    }

    /**
       Returns the arc tangent of an angle. The results may have
       maximum error of 2 ulps.
     */
    public static double atan(double s) {
        double t, u;
        int q = 0;

        if (s < 0) { s = -s; q = 2; }
        if (s > 1) { s = 1.0 / s; q |= 1; }

        t = s * s;

        u = -1.88796008463073496563746e-05;
        u = u * t + (0.000209850076645816976906797);
        u = u * t + (-0.00110611831486672482563471);
        u = u * t + (0.00370026744188713119232403);
        u = u * t + (-0.00889896195887655491740809);
        u = u * t + (0.016599329773529201970117);
        u = u * t + (-0.0254517624932312641616861);
        u = u * t + (0.0337852580001353069993897);
        u = u * t + (-0.0407629191276836500001934);
        u = u * t + (0.0466667150077840625632675);
        u = u * t + (-0.0523674852303482457616113);
        u = u * t + (0.0587666392926673580854313);
        u = u * t + (-0.0666573579361080525984562);
        u = u * t + (0.0769219538311769618355029);
        u = u * t + (-0.090908995008245008229153);
        u = u * t + (0.111111105648261418443745);
        u = u * t + (-0.14285714266771329383765);
        u = u * t + (0.199999999996591265594148);
        u = u * t + (-0.333333333333311110369124);

        t = s + s * (t * u);

        if ((q & 1) != 0) t = 1.570796326794896557998982 - t;
        if ((q & 2) != 0) t = -t;

        return t;
    }

    private static final double PI4_A = .7853981554508209228515625;
    private static final double PI4_B = .794662735614792836713604629039764404296875e-8;
    private static final double PI4_C = .306161699786838294306516483068750264552437361480769e-16;
    private static final double M_1_PI = 0.3183098861837906715377675267450287;

    /**
       Returns the trigonometric sine of an angle. The results may
       have maximum error of 2 ulps.
     */
    public static double sin(double d) {
        int q;
        double u, s;

        u = d * M_1_PI;
        q = (int)(u < 0 ? u - 0.5 : u + 0.5);

        d = mla(q, -PI4_A*4, d);
        d = mla(q, -PI4_B*4, d);
        d = mla(q, -PI4_C*4, d);

        if ((q & 1) != 0) d = -d;

        s = d * d;

        u = -7.97255955009037868891952e-18;
        u = mla(u, s, 2.81009972710863200091251e-15);
        u = mla(u, s, -7.64712219118158833288484e-13);
        u = mla(u, s, 1.60590430605664501629054e-10);
        u = mla(u, s, -2.50521083763502045810755e-08);
        u = mla(u, s, 2.75573192239198747630416e-06);
        u = mla(u, s, -0.000198412698412696162806809);
        u = mla(u, s, 0.00833333333333332974823815);
        u = mla(u, s, -0.166666666666666657414808);

        u = mla(s, u * d, d);

        return u;
    }

    /**
       Returns the trigonometric cosine of an angle. The results may
       have maximum error of 2 ulps.
     */
    public static double cos(double d) {
        int q;
        double u, s;

        q = 1 + 2*rint(d * M_1_PI - 0.5);

        d = mla(q, -PI4_A*2, d);
        d = mla(q, -PI4_B*2, d);
        d = mla(q, -PI4_C*2, d);

        if ((q & 2) == 0) d = -d;

        s = d * d;

        u = -7.97255955009037868891952e-18;
        u = mla(u, s, 2.81009972710863200091251e-15);
        u = mla(u, s, -7.64712219118158833288484e-13);
        u = mla(u, s, 1.60590430605664501629054e-10);
        u = mla(u, s, -2.50521083763502045810755e-08);
        u = mla(u, s, 2.75573192239198747630416e-06);
        u = mla(u, s, -0.000198412698412696162806809);
        u = mla(u, s, 0.00833333333333332974823815);
        u = mla(u, s, -0.166666666666666657414808);

        u = mla(s, u * d, d);

        return u;
    }
}
