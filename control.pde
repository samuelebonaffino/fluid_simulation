import controlP5.*;

ControlP5 cp5_dMul, cp5_vMul;
int startY = 10, incrY = 30;

void initControlSystem()
{
    cp5_dMul = initControl("dMul", startY, 1, 150, 0, 300);
    cp5_vMul = initControl("vMul", startY+incrY, 1, 2, -150, 150);
}

ControlP5 initControl(String name, int y, float mul, float start)
{
    ControlP5 cp5 = new ControlP5(this);
    cp5.addNumberbox(name)
       .setPosition(10, y)
       .setSize(30, 15)
       .setRange(-1, 1)
       .setMultiplier(mul)
       .setDirection(Controller.HORIZONTAL)
       .setValue(start);
    return cp5;
}
ControlP5 initControl(String name, int y, float mul, float start, int min, int max)
{
    ControlP5 cp5 = new ControlP5(this);
    cp5.addNumberbox(name)
       .setPosition(10, y)
       .setSize(30, 15)
       .setRange(min, max)
       .setMultiplier(mul)
       .setDirection(Controller.HORIZONTAL)
       .setValue(start);
    return cp5;
}