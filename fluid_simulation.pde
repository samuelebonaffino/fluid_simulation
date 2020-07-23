import processing.sound.*;

final static int N = 128;
final static int SCALE = 4;
final static int I = 16;

Fluid fluid;
Audio audio;
float t = 0;
float dMul, vMul;

void settings()
{
    size(N*SCALE, N*SCALE);
}

void setup() 
{
    initControlSystem();
    fluid = new Fluid(0.2, 0, 0.0000001);  
    audio = new Audio(128, "bias.wav");
    audio.play();
}

void draw() 
{
    background(0);

    addDensity(audio);
    addVelocity(audio);

    fluid.step(audio);
    fluid.renderD();
    fluid.fadeD();
}

void addDensity()
{
    int cx = int(0.5*width/SCALE);
    int cy = int(0.5*height/SCALE);
    for(int i = -1; i <= 1; i++)
        for(int j = -1; j <= 1; j++)
            fluid.addDensity(cx+i, cy+j, random(50, 100));
}
void addDensity(Audio audio)
{
    int cx = int(0.5*width/SCALE);
    int cy = int(0.5*height/SCALE);
    for(int i = -1; i <= 1; i++)
        for(int j = -1; j <= 1; j++)
            fluid.addDensity(cx+i, cy+j, audio.getAmplitude(dMul));
}

void addVelocity()
{
    int cx = int(0.5*width/SCALE);
    int cy = int(0.5*height/SCALE);
    for(int i = 0; i < 2; i++)
    {
        float angle = noise(t)*2*TWO_PI;
        PVector v = PVector.fromAngle(angle);
        v.mult(0.2);
        t += 0.01;
        fluid.addVelocity(cx, cy, v.x, v.y);
    }
}
void addVelocity(Audio audio)
{
    int cx = int(0.5*width/SCALE);
    int cy = int(0.5*height/SCALE);
    for(int i = 0; i < 2; i++)
    {
        float angle = noise(t)*2*TWO_PI;
        PVector v = PVector.fromAngle(angle);
        v.mult(audio.getAmplitude(vMul));
        t += 0.01;
        fluid.addVelocity(cx, cy, v.x, v.y);
    }
}