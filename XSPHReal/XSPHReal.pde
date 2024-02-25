
String projectTitle = "SPH Fluid Sim";


float floor = 695;
float ceiling = 0;
float left = 0;
float right = 1280;

float delta = 20;
float dtime = .01;
//float gravity = 9.8;
//float pull = .1;

float pointsize = 22;

float smooth = 17;
float stiff = 50000;
float rest = 5.0;

int numpoints = 4000;
float grablength = 50;

class Particle
{
  PVector pos;
  PVector oldPos;
  PVector vel;
  float dens, press;
  int id;
  
  Particle(float x, float y, int ID)
  {
    oldPos = new PVector(x,y);
    pos = new PVector(x,y);
    vel = new PVector();
    dens = 0;
    press = 0;
    id = ID;
  }
  
  void drawParticle(){
    float pressred = 30;
    float pressblue = 2;
    color dye;
    float diff = pressred - pressblue;
    float quad = diff / 4;
    if (press < pressblue) dye = color(50, 50, 255); 
    else if(press < pressblue+quad){
      float segdiff = press-pressblue;
      dye = color(50, 255*(segdiff/quad), 255);
    }
    else if(press < pressblue+2*quad){
      float segdiff = press-(pressblue+quad);
      dye = color(50, 255, 255-255*(segdiff/quad));
    }
    else if(press < pressblue+3*quad){
      float segdiff = press-(pressblue+2*quad);
      dye = color(255*(segdiff/quad), 255, 50);
    }
    else if(press < pressred){
      float segdiff = press-(pressblue+3*quad);
      dye = color(255, 255-255*(segdiff/quad), 50);
    }
    else dye = color(255, 50, 50);
    stroke(dye);
    point(pos.x, pos.y);
  }
  
}



class Pair{
  Particle a, b;
  float x, x2, x3;
  Pair(Particle A, Particle B, float dist) {
    a = A; b = B;
    x = dist/smooth;
    x2 = pow(x,2);
    x3 = pow(x,3);
  }
}


class SphFluid 
{

  ArrayList<Particle> particles;
  ArrayList<Pair> pairs;
  int num;
  float ksmooth, kstiff, kstiffN, krest;
  float reach;
  float WConstant;
  float gradWConstant;
  float alpha;
  
  PVector oldMouse;
  PVector Mouse;
  
  SphFluid(int num, float ksm, float kst, float kr, float reach) {
    this.num = num;
    ksmooth = ksm;
    kstiff = kst;
    krest = kr;
    this.reach = reach;
    particles = new ArrayList<Particle>();
    WConstant = PI*pow(ksmooth,3);
    gradWConstant = PI*pow(ksmooth,4);
    alpha = 0.3;
    oldMouse = new PVector(mouseX, mouseY);
    Mouse = new PVector(mouseX, mouseY);
    
    float interval = 15.0*0.97;
    float initx = interval;
    float inity = interval;
    
    for (int i = 0; i < num; i++) {
      particles.add(new Particle(initx, inity, i));
      initx += interval;
      if (initx > right - interval) {
        initx = interval;
        inity += interval;
      }
    }
    
    pairs = new ArrayList<Pair>();
  }
  
  PVector displacement(float x, float y, float U, float L){
    float r = dist(0,0,x,y);
    float eps = 0.00001;
    PVector v = new PVector(r*L-y*y, x*y);
    v.mult(U/(r*L*exp(r/(L+eps))+eps));
    return v;
  }
  
  void updateParticles(float dt){
      pairs.clear();
      
      Mouse.set(mouseX, mouseY);
      PVector ctrl = PVector.sub(Mouse,oldMouse);
      PVector midMouse = PVector.add(oldMouse, PVector.mult(ctrl, 0.5));
      float lambda = constrain(ctrl.mag(),0,100);
      ctrl.normalize();
      PVector spinx = new PVector(ctrl.x, -ctrl.y);
      PVector spiny = new PVector(ctrl.y, ctrl.x);
      PVector diff = new PVector();
      
      for(int i=0; i<num; i++)
      {
          Particle p = particles.get(i);
          p.vel = PVector.sub(p.pos, p.oldPos).div(dt); // get old velocity
        
          if(mousePressed && mouseButton==RIGHT){
              diff = PVector.sub(p.pos, oldMouse);
              float x = ctrl.dot(diff);
              float y = ctrl.cross(diff).z;
              PVector v = displacement(x, y, 0.2*lambda/dt, 15).mult(0.5);
              p.vel.x += spinx.dot(v);
              p.vel.y += spiny.dot(v);
          }
        
          if(mousePressed && mouseButton==RIGHT){
              diff = PVector.sub(p.pos, midMouse);
              float x = ctrl.dot(diff);
              float y = ctrl.cross(diff).z;
              PVector v = displacement(x, y, 0.2*lambda/dt, 15).mult(0.5);
              p.vel.x += spinx.dot(v);
              p.vel.y += spiny.dot(v);
          }
        
          if(p.pos.x<=left+pointsize/2){ p.pos.x=left+(pointsize/2)*1.05; p.vel.x *= -0.3; }
          else if(p.pos.x>=(right-pointsize/2)){ p.pos.x = right-(pointsize/2)*1.05; p.vel.x *= -0.3; }
          if(p.pos.y>=floor-pointsize/2){ p.pos.y=floor-pointsize/2*1.05; p.vel.y *= -0.2; }
          else if(p.pos.y<=ceiling+pointsize/2){ p.pos.y=ceiling+pointsize/2*1.05; p.vel.y *= -0.2; }
          p.oldPos.set(p.pos);
          p.pos.add(p.vel.copy().mult(dt));
          p.dens = 0; p.press = 0;
      }
      oldMouse.set(Mouse);
      
      for(int i=0; i<num; i++){
        for(int j=0; j<num; j++){
          if(i != j){
            Particle p1 = particles.get(i);
            Particle p2 = particles.get(j);
            float dist = dist(p1.pos.x, p1.pos.y, p2.pos.x, p2.pos.y); 
            if(dist < 2*ksmooth) pairs.add(new Pair(p1, p2, dist));
          }
        }
      }
      
      for(int i=0; i<pairs.size(); i++){
        Pair p = pairs.get(i);
        float density = 1.0*W(p.x, p.x2, p.x3);
        p.a.dens += density;
        p.b.dens += density;
      }
      
      for(int i=0; i<num; i++){
        Particle p = particles.get(i);
        p.press = kstiff*p.dens-krest;
        p.press = constrain(p.press, 0, 30);
      }
      
      for(int i=0; i<pairs.size(); i++){
        Pair p = pairs.get(i);
        float press = p.a.press/pow(p.a.dens,2) + p.b.press/pow(p.b.dens,2);
        float displace = -(press*gradW(p.x, p.x2))*dt;
        PVector ab = PVector.sub(p.a.pos, p.b.pos).normalize().mult(displace*dt);
        p.a.pos.add(ab);
        p.b.pos.sub(ab);
      }
      
      for(int i = 0; i < num; i++){
         Particle p = particles.get(i);
         p.vel = PVector.sub(p.pos, p.oldPos).div(dt);
      }
      
      for(int i=0; i<pairs.size(); i++){
        Pair p = pairs.get(i);
        float constant = alpha*0.5*(1.0/p.a.dens + 1.0/p.b.dens);
        PVector ab = PVector.sub(p.b.vel, p.a.vel).mult(W(p.x, p.x2, p.x3)*constant*dt);
        p.a.pos.add(ab);
        p.b.pos.sub(ab);
      }
    
  }
  
  
  float W(float x, float x2, float x3){ // x = dist/ksmooth
     return ((x<1.0) ? (1.0-1.5*x2+0.75*x3) : (0.25*pow(2-x,3))) / (WConstant);
  }
  
  
  float gradW(float x, float x2){
     return ((x<1.0) ? (2.25*x2-3.0*x) : (-0.75*pow(2-x,2))) / (gradWConstant);
  }
  
  void drawParticles() {
    strokeWeight(pointsize);
    for (int i = 0; i < num; i++) {
      Particle p = particles.get(i);
      p.drawParticle();
    }
    strokeWeight(2); stroke(0);
    for(int i=0; i<pairs.size(); i++){
      Pair p = pairs.get(i);
      if(p.a==particles.get(1500)) line(p.a.pos.x, p.a.pos.y, p.b.pos.x, p.b.pos.y);
    }
    noStroke();
  }
}




SphFluid fluid = new SphFluid(numpoints, smooth, stiff, rest, grablength);


void setup() {
 size(1280, 695, P3D);
 noStroke();
}

void computePhysics(float dt) {
  fluid.updateParticles(dt);
}

void drawScene(){
  background(50, 51, 54);
  fluid.drawParticles();
}

void draw() {
  float startFrame = millis(); 
  computePhysics(dtime); 
  float endPhysics = millis();
  
  drawScene();
  float endFrame = millis();
  delta = endFrame - startFrame;
  
  String runtimeReport = "Frame: "+str(endFrame-startFrame)+"ms,"+
        " Physics: "+ str(endPhysics-startFrame)+"ms,"+
        " FPS: "+ str(round(frameRate)) + "\n";
  surface.setTitle(projectTitle+ "  -  " +runtimeReport);
}
