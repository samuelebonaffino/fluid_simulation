final int N = 128;
final int SCALE = 4;
final int I = 16;

Fluid fluid;
float t = 0;

void settings()
{
    size(N*SCALE, N*SCALE);
}

void setup() 
{
    fluid = new Fluid(0.2, 0, 0.0000001);    
}

void draw() 
{
    background(0);

    int cx = int(0.5*width/SCALE);
    int cy = int(0.5*height/SCALE);
    for(int i = -1; i <= 1; i++)
        for(int j = -1; j <= 1; j++)
            fluid.addDensity(cx+i, cy+j, random(50, 100));

    for(int i = 0; i < 2; i++)
    {
        float angle = noise(t)*2*TWO_PI;
        PVector v = PVector.fromAngle(angle);
        v.mult(0.2);
        t += 0.01;
        fluid.addVelocity(cx, cy, v.x, v.y);
    }

    fluid.step();
    fluid.renderD();
    fluid.fadeD();
}
