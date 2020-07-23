class Fluid
{   
    float dt, diff, visc;
    float[] s;
    float[] density;
    float[] vx;
    float[] vy;
    float[] vx0;
    float[] vy0;

    Fluid(float dt, float diff, float visc)
    {
        this.dt   = dt;
        this.diff = diff;
        this.visc = visc;

        s         = new float[N*N];
        density   = new float[N*N];
        vx        = new float[N*N];
        vy        = new float[N*N];
        vx0       = new float[N*N];
        vy0       = new float[N*N];
    }

    void addDensity(int x, int y, float amount)
    {
        density[IX(x,y)] += amount;
    }

    void addVelocity(int x, int y, float amountX, float amountY)
    {
        vx[IX(x,y)] += amountX;
        vy[IX(x,y)] += amountY;
    }

    void step()
    {
        diffuse(1, vx0, vx, visc, dt);
        diffuse(2, vy0, vy, visc, dt);

        project(vx0, vy0, vx, vy);

        advect(1, vx, vx0, vx0, vy0, dt);
        advect(2, vy, vy0, vx0, vy0, dt);
        
        project(vx, vy, vx0, vy0);

        diffuse(0, s, density, diff, dt);
        advect(0, density, s, vx, vy, dt);
    }

    void renderD()
    {
        colorMode(HSB, 255);

        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
            {
                float x = i * SCALE;
                float y = j * SCALE;
                float d = this.density[IX(i,j)];
                fill((d+50)%255, 200, d);
                noStroke();
                rect(x, y, SCALE, SCALE);
            }
    }

    void renderV()
    {
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
            {
                float x = i * SCALE;
                float y = j * SCALE;
                float vx = this.vx[IX(i,j)];
                float vy = this.vy[IX(i,j)];
                stroke(255);

                if(!(abs(vx) < 0.1 && abs(vy) <= 0.1))
                    line(x, y, x + vx*SCALE, y + vy*SCALE);
            }
    }

    void fadeD()
    {
        for(int i = 0; i < this.density.length; i++)
        {
            float d = this.density[i];
            this.density[i] = constrain(d - 0.02, 0, 255);
            // this.density[i] *= 0.99;
        }
    }
}