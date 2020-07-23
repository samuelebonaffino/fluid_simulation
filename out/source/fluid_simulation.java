import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class fluid_simulation extends PApplet {

final int N = 128;
final int SCALE = 4;
final int I = 16;

Fluid fluid;
float t = 0;

public void settings()
{
    size(N*SCALE, N*SCALE);
}

public void setup() 
{
    fluid = new Fluid(0.2f, 0, 0.0000001f);    
}

public void draw() 
{
    background(0);

    int cx = PApplet.parseInt(0.5f*width/SCALE);
    int cy = PApplet.parseInt(0.5f*height/SCALE);
    for(int i = -1; i <= 1; i++)
        for(int j = -1; j <= 1; j++)
            fluid.addDensity(cx+i, cy+j, random(50, 100));

    for(int i = 0; i < 2; i++)
    {
        float angle = noise(t)*2*TWO_PI;
        PVector v = PVector.fromAngle(angle);
        v.mult(0.2f);
        t += 0.01f;
        fluid.addVelocity(cx, cy, v.x, v.y);
    }

    fluid.step();
    fluid.renderD();
    fluid.fadeD();
}
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

    public void addDensity(int x, int y, float amount)
    {
        density[IX(x,y)] += amount;
    }

    public void addVelocity(int x, int y, float amountX, float amountY)
    {
        vx[IX(x,y)] += amountX;
        vy[IX(x,y)] += amountY;
    }

    public void step()
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

    public void renderD()
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

    public void renderV()
    {
        for(int i = 0; i < N; i++)
            for(int j = 0; j < N; j++)
            {
                float x = i * SCALE;
                float y = j * SCALE;
                float vx = this.vx[IX(i,j)];
                float vy = this.vy[IX(i,j)];
                stroke(255);

                if(!(abs(vx) < 0.1f && abs(vy) <= 0.1f))
                    line(x, y, x + vx*SCALE, y + vy*SCALE);
            }
    }

    public void fadeD()
    {
        for(int i = 0; i < this.density.length; i++)
        {
            // float d = this.density[i];
            // this.density[i] = constrain(d - 0.02, 0, 255);
            this.density[i] *= 0.99f;
        }
    }
}
public int IX(int x, int y)
{
    x = constrain(x, 0, N-1);
    y = constrain(y, 0, N-1);
    return x + y*N;
}

public void diffuse(int b, float[] x, float[] x0, float diff, float dt)
{
    float a = dt*diff*(N-2)*(N-2);
    lin_solve(b, x, x0, a, 1+4*a);
}

public void lin_solve(int b, float[] x, float[] x0, float a, float c)
{
    float cRecip = 1.0f/c;
    for(int k = 0; k < I; k++)
        for(int j = 1; j < N-1; j++)
            for(int i = 1; i < N-1; i++)
            {
                x[IX(i,j)] = (x0[IX(i,j)]
                              + a*(x[IX(i+1,j)]
                              + x[IX(i-1,j)]
                              + x[IX(i,j+1)]
                              + x[IX(i,j-1)]
                              )) * cRecip;
            }
    set_bnd(b, x);
}

public void project(float[] velocX, float[] velocY, float[] p, float[] div)
{
    for(int j = 1; j < N-1; j++)
        for(int i = 1; i < N-1; i++)
        {
            div[IX(i,j)] = -0.5f*(
                  velocX[IX(i+1,j)]
                - velocX[IX(i-1,j)]
                + velocY[IX(i,j+1)]
                - velocY[IX(i,j-1)]
                )/N;
                p[IX(i,j)] = 0;
        }
    
    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 4);

    for(int j = 1; j < N-1; j++)
        for(int i = 1; i < N-1; i++)
        {
            velocX[IX(i,j)] -= 0.5f*(p[IX(i+1,j)]-p[IX(i-1,j)])*N;
            velocY[IX(i,j)] -= 0.5f*(p[IX(i,j+1)]-p[IX(i,j-1)])*N;
        }
    
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}

public void advect(int b, float[] d, float[] d0, float[] velocX, float[] velocY, float dt)
{
    float i0, i1, j0, j1;
    float dtx = dt*(N-2);
    float dty = dt*(N-2);
    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;
    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;

    for(j = 1, jfloat = 1; j < N-1; j++, jfloat++)
        for(i = 1, ifloat = 1; i < N-1; i++, ifloat++)
        {
            tmp1 = dtx*velocX[IX(i,j)];
            tmp2 = dty*velocY[IX(i,j)];
            x = ifloat-tmp1;
            y = jfloat-tmp2;

            if(x < 0.5f)         x = 0.5f;
            if(x > Nfloat+0.5f)  x = Nfloat+0.5f;
            i0 = floor(x);
            i1 = i0+1.0f;
            if(y < 0.5f)         y = 0.5f;
            if(y > Nfloat+0.5f)  y = Nfloat+0.5f;
            j0 = floor(y);
            j1 = j0+1.0f;

            s1 = x-i0;
            s0 = 1.0f-s1;
            t1 = y-j0;
            t0 = 1.0f-t1;

            int i0i = PApplet.parseInt(i0);
            int i1i = PApplet.parseInt(i1);
            int j0i = PApplet.parseInt(j0);
            int j1i = PApplet.parseInt(j1);

            d[IX(i,j)] = s0*(t0*d0[IX(i0i,j0i)]+t1*d0[IX(i0i,j1i)]) +
                         s1*(t0*d0[IX(i1i,j0i)]+t1*d0[IX(i1i,j1i)]);
                         
        }
        set_bnd(b, d);
}

public void set_bnd(int b, float[] x)
{
    for(int i = 1; i < N-1; i++)
    {
        x[IX(i,0)]   = b == 2 ? -x[IX(i,1)] : x[IX(i,1)];
        x[IX(i,N-1)] = b == 2 ? -x[IX(i,N-2)] : x[IX(i,N-2)];
    }
    for(int j = 1; j < N-1; j++)
    {
        x[IX(0,j)]   = b == 1 ? -x[IX(1,j)] : x[IX(1,j)];
        x[IX(N-1,j)] = b == 1 ? -x[IX(N-2,j)] : x[IX(N-2,j)];
    }

    x[IX(0,0)]     = 0.5f*(x[IX(1,0)])+x[IX(0,1)];
    x[IX(0,N-1)]   = 0.5f*(x[IX(1,N-1)])+x[IX(0,N-2)];
    x[IX(N-1,0)]   = 0.5f*(x[IX(N-2,0)])+x[IX(N-1,1)];
    x[IX(N-1,N-1)] = 0.5f*(x[IX(N-2,N-1)])+x[IX(N-1,N-2)];
}
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "fluid_simulation" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}