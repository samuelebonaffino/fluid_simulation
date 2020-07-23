class Fluid()
{
    final static int SIZE = 8;
    
    float diff, visc;
    float[] s;
    float[] density;
    float[] vx;
    float[] vy;
    float[] vx0;
    float[] vy0;

    Fluid(float diff, float visc)
    {
        int N = SIZE;

        this.diff = diff;
        this.visc = visc;

        s         = new float[N*N];
        density   = new float[N*N];
        vx        = new float[N*N];
        vy        = new float[N*N];
        vx0       = new float[N*N];
        vy0       = new float[N*N];
    }

    int getLinearIndex(int x, int y)
    {
        return x + y*SIZE;
    }

    void addDensity(int x, int y, float amount)
    {
        density[getLinearIndex(x,y)] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {
        vx[getLinearIndex(x,y)] += amountX;
        vy[getLinearIndex(x,y)] += amountY;
    }

    
}