
/*********************************************************/
static float smooth = 17; 
static float dt = 0.025;
static float size = 22;
static float stiff = 55000;
static float P0 = 5.0;
static float alpha = 0.3;
static float Wconst = 1.0/(PI*pow(smooth,3));
static float gradWconst = 1.0/(PI*pow(smooth,4));

color Rainbow(float x){
    float RED = 30;
    float BLUE = 2;
    float q = (RED-BLUE)/4.0;
    if(x < BLUE) return color(50,50,255);
    else if(x < BLUE+q) return color(50, 255*(x-BLUE)/q, 255);
    else if(x < BLUE+2*q) return color(50, 255, 255*2-255*(x-BLUE)/q);
    else if(x < RED-q) return color(255*(x-BLUE)/q-255*2, 255, 50);
    else if(x < RED) return color(255, 255*4-255*(x-BLUE)/q, 50);
    return color(255, 50, 50);
}
/*********************************************************/

class Particle{
  PVector pos, vel, acc;
  float density;
  float pressure;
  int id;
  
  Particle(float x, float y, int id_){
    id = id_;
    pos = new PVector(x,y);
    vel = new PVector();
    acc = new PVector();
    density = 0;
    pressure = 0;
  }
  
  void display(){
    stroke(Rainbow(pressure));
    stroke(Rainbow((float)id/140));
    point(pos.x, pos.y);
  }
}


class Pair{
  Particle a, b;
  float x, x2, x3;
  Pair(Particle a_, Particle b_, float dist){
    a = a_; b = b_;
    x = dist/smooth;
    x2 = pow(x,2);
    x3 = pow(x,3);
  }
}

class Mouse{
  PVector oldPos, pos;
  Mouse(){
    oldPos = new PVector();
    pos = new PVector();
  }
  void update(){
    oldPos.set(pos);
    pos.set(mouseX, mouseY);
  } 
}

/***********************************************************/

class Fluid{
  
  ArrayList<Particle> particles;
  ArrayList<Pair> pairs; 
  int N;
  int debt;
  Mouse mouse;
  Palette board;
  
  Fluid(int n_){
    N = n_; debt = 0;
    particles = new ArrayList<Particle>();
    pairs = new ArrayList<Pair>();
    mouse = new Mouse();
    board = new Palette(this);
    createFluid();
  }
  
  void createFluid(){
    float interval = 15*0.97;
    int x=1, y=1;
    for(int i=0; i<N; i++){
      particles.add(new Particle(interval*x+randomGaussian(), interval*y+randomGaussian(), i));
      x++; if(interval*x > width){ x=1; y++; } 
    }
  }
  
  void update(){
    pairs.clear();
    mouse.update();
    if(mousePressed && mouseButton==LEFT) mouseDrag(); 
    for(Particle p : particles){ p.density = 0; p.pressure = 0; p.acc = new PVector(); }
    addPairs();
    calculate();
    SPH();  board.update();  moveParticles();
    XSPH(); moveParticles();
    checkEdge();
  }
  
  void mouseDrag(){ // add to pos : no acceleration, just deformation effect
    int steps = 5;
    PVector pos = mouse.oldPos.copy();
    PVector diff = PVector.sub(mouse.pos, mouse.oldPos);
    float lambda = constrain(diff.mag(), 0, 100);
    PVector dir = diff.copy().normalize();
    PVector spinx = new PVector(dir.x, -dir.y);
    PVector spiny = new PVector(dir.y, dir.x);
    for(int j=0; j<steps; j++){
      pos.add(PVector.div(diff,(float)steps));
      for(Particle p : particles){
        PVector ctrl = PVector.sub(p.pos, pos);
        float x = dir.dot(ctrl);
        float y = dir.cross(ctrl).z;
        PVector V = displacement(x, y, 0.1*lambda/dt, 20).div((float)steps);
        PVector D = displacement(x, y, lambda/dt, 25).div((float)steps).mult(dt);
        p.vel.x += spinx.dot(V);
        p.vel.y += spiny.dot(V);
      }
      for(ArrayList<System> line : board.lines){
        for(System s : line){
          ColorPoint p = s.point;
          PVector ctrl = PVector.sub(p.pos, pos);
          float x = dir.dot(ctrl);
          float y = dir.cross(ctrl).z;
          PVector V = displacement(x, y, 0.1*lambda/dt, 20).div((float)steps);
          PVector D = displacement(x, y, lambda/dt, 25).div((float)steps).mult(dt);
          p.vel.x += spinx.dot(V);
          p.vel.y += spiny.dot(V);
          p.pos.x += spinx.dot(D);
          p.pos.y += spiny.dot(D);
        }
      }
    }
  }
  
  void addDrop(float x, float y){
    int nums = 100;
    for(int i=0; i<nums; i++){ particles.add(new Particle(x-20+random(40),y-20+random(40),N+debt)); debt++; }
    board.addDrop(x,y);
  }
  
  void checkEdge(){
    for(int i=N+debt-1; i>=0; i--){
      Particle p = particles.get(i);
      if(p.pos.x<=size/2 || p.pos.x>=width-size/2){
        if(debt>0){ debt--; particles.remove(i); continue; }
        p.vel.x *= -0.3;
        if(p.pos.x<=size/2) p.pos.x = (size/2)*1.05;
        else p.pos.x = width-(size/2)*1.05;
      }
      if(p.pos.y<=size/2 || p.pos.y>=height-size/2){
        if(debt>0){ debt--; particles.remove(i); continue; }
        p.vel.y *= -0.3;
        if(p.pos.y<=size/2) p.pos.y = (size/2)*1.05;
        else p.pos.y = height-(size/2)*1.05;
      }
    }
  }
  
  void addPairs(){
    for(int i=0; i<N; i++){
      Particle p1 = particles.get(i);
      for(int j=0; j<N; j++) if(i!=j){
        Particle p2 = particles.get(j);
        float dist = dist(p1.pos.x, p1.pos.y, p2.pos.x, p2.pos.y);
        if(dist<2*smooth) pairs.add(new Pair(p1, p2, dist));
      }
    }
  }
  
  void calculate(){
    for(Pair p : pairs){
      float d = W(p.x, p.x2, p.x3); // let's say m = 1.0
      p.a.density += d;
      p.b.density += d;
    }
    for(Particle p : particles)
      p.pressure = constrain(stiff*p.density-P0, 0, 30);
  }
  
  void SPH(){
    for(Pair p : pairs){
      float P = p.a.pressure/pow(p.a.density,2) + p.b.pressure/pow(p.b.density,2);
      float acc = -(P*gradW(p.x, p.x2));
      PVector dir = PVector.sub(p.a.pos, p.b.pos).normalize();
      PVector dv = PVector.mult(dir,acc*dt);
      p.a.vel.add(dv);
      p.b.vel.sub(dv);
    }
  }
  
  void XSPH(){
    for(Pair p : pairs){
      float constant = alpha*0.5*(1.0/p.a.density + 1.0/p.b.density);
      PVector dv = PVector.sub(p.b.vel, p.a.vel).mult(W(p.x, p.x2, p.x3)*constant);
      p.a.vel.add(dv);
      p.b.vel.sub(dv);
    }
  }
  
  void moveParticles(){
    for(Particle p : particles){
      p.pos.add(p.vel.copy().mult(dt));
    }
  }
  
  void display(){
    strokeWeight(size);
    if(keyPressed && key=='m') for(Particle p : particles) p.display();
    board.display();
  }
  
  PVector displacement(float x, float y, float U, float L){
    float r = dist(0,0,x,y);
    float eps = 0.00001;
    PVector v = new PVector(r*L-y*y, x*y);
    v.mult(U/(r*L*exp(r/(L+eps))+eps));
    return v;
  }
  
  float W(float x, float x2, float x3){
    return ((x<1.0) ? (1.0-1.5*x2+0.75*x3):(0.25*pow(2-x,3))) * Wconst;
  }
  
  float gradW(float x, float x2){
    return ((x<1.0) ? (2.25*x2-3.0*x):(-0.75*pow(2-x,2))) * gradWconst;
  }
  
}

/*************************************************************/

static float dx=1;
static float SCALE=5;
static float vSmooth = 17;
static float vWconst = 1.0/(PI*pow(vSmooth,3));

class Palette{
  ArrayList<ArrayList<System>> lines;
  ArrayList<ArrayList<Integer>> cut;
  Fluid fluid;
  //Boolean[] ring;
  color[] colors;
  
  Palette(Fluid fluid_){
    lines = new ArrayList();
    cut = new ArrayList();
    colors = new color[0];
    fluid = fluid_;
  }
  
  void update(){
    addPairs();
    calculate();
    movePoints();
    addRemove();
  }
  
  void addDrop(float X, float Y){
    int n = 30;
    float r = 50;
    color dye = color(100+random(100), 100+random(100), 100+random(100));
    ArrayList<System> line = new ArrayList();
    for(int i=0; i<=n; i++) line.add(new System(X+r*cos(2*PI*i/n+1),Y+r*sin(2*PI*i/n+1)));
    lines.add(line); colors = append(colors, dye);
  }
  
  void addPairs(){
    for(ArrayList<System> line : lines)
      for(System s1 : line){
        ColorPoint p1 = s1.point;
        p1.pair.clear();
        for(Particle p2 : fluid.particles){
          float dist = dist(p1.pos.x, p1.pos.y, p2.pos.x, p2.pos.y);
          if(dist<2*vSmooth) p1.pair.add(p2);
        }
      }
  }
  
  void calculate(){
    for(ArrayList<System> line : lines)
      for(System s1 : line){
        ColorPoint p1 = s1.point;
        PVector sum = new PVector();
        int count=0;
        //PVector max = new PVector();
        for(Particle p2 : p1.pair){ 
          PVector vel = p2.vel.copy();
          //sum.add(vel); count++; 
          //if(vel.mag()>max.mag()) max.set(vel);
          float dist = dist(p1.pos.x, p1.pos.y, p2.pos.x, p2.pos.y);
          sum.add(vel.mult(W(dist/vSmooth)/(p2.density+0.00001))); 
          //sum.add(p2.vel.copy().div(p1.pair.size()));
        }
        if(count!=0) sum.div(count);
        p1.vel.add(sum.mult(1.0));
      }
  }
  
  void movePoints(){
    for(ArrayList<System> line : lines)
      for(System s : line){
        ColorPoint p = s.point;
        p.update();
        p.vel = new PVector();
      }
  }
  
  void display(){
    
    for(int I=0; I<lines.size(); I++){
      ArrayList<System> line = lines.get(I);
      fill(colors[I]);
      int n = line.size()-1;
      line.get(0).setValue(line.get(1));
      for(int i=1; i<n-1; i++)
        line.get(i).setValue(line.get(i-1), line.get(i+1));
      line.get(n-1).setValue(line.get(n-2), line.get(n), 0);
      line.get(n-1).solve(line.get(n), 0);
      for(int i=n-2; i>=0; i--)
        line.get(i).solve(line.get(i+1));
      
      noStroke();
      beginShape();
      for(int i=0; i<n; i++)
        line.get(i).show(line.get(i+1));
      endShape();
      
      stroke(150);
      strokeWeight(3);
      if(keyPressed && key==' ') for(int i=0; i<line.size()-1; i++){
        System s1 = line.get(i);
        if(i==0) stroke(150, 0, 0); 
        else stroke(250);
        s1.point.display(); 
      }
      
    }
  }
  
  void addRemove(){
    for(ArrayList<System> line : lines){
      int n = line.size()-1;
      for(int i=0; i<n; i++){
        System s1 = line.get(i);
        System s2 = line.get(i+1);
        if(PVector.sub(s1.point.pos, s2.point.pos).mag()<5){
          line.remove(i+1); i--; n--;
        }
      }
      for(int i=0; i<n; i++){
        System s1 = line.get(i);
        System s2 = line.get(i+1);
        if(PVector.sub(s1.point.pos, s2.point.pos).mag()>30){
          System m = s1.midPoint(s2);
          line.add(i+1, m); i++; n++;
        }
      }
    }
  }
  
  float W(float x){
    return ((x<1.0) ? (1.0-1.5*pow(x,2)+0.75*pow(x,3)):(0.25*pow(2-x,3))) * vWconst;
  }
  
}

class ColorPoint{
   PVector pos, vel;
   ArrayList<Particle> pair;
   ColorPoint(float X, float Y){
      pos = new PVector(X,Y);
      vel = new PVector();
      pair = new ArrayList<Particle>();
   }
   void update(){ pos.add(vel.copy().mult(dt)); }
   void display(){ point(pos.x, pos.y); }
}

class System{
  PVector a, b; // need to solve
  PVector r, R; 
  float Z; 
  ColorPoint point;
  static final int n = 10;
  
  System(float X, float Y){
    a = new PVector();
    b = new PVector();
    r = new PVector();
    R = new PVector();
    point = new ColorPoint(X,Y);
  }
  
  PVector getPos(){ return point.pos.copy(); }
  
  void setValue(System next){
    r = getPos().add(next.getPos().mult(2.0)); 
    R = PVector.div(r,2.0);
    Z = 0.5;
  }
  void setValue(System prev, System next){
    r = getPos().mult(2.0).add(next.getPos()).mult(2.0);
    R = PVector.sub(r,prev.R).div(4.0-prev.Z);
    Z = 1/(4.0-prev.Z);
  }
  void setValue(System prev, System next, int last){
    r = getPos().mult(8.0).add(next.getPos());
    R = PVector.sub(r, PVector.mult(prev.R,2.0)).div(7.0-2.0*prev.Z);
  }
  
  void solve(System next, int last){
    a = R.copy();
    b = PVector.add(a, next.getPos()).div(2.0);
  }
  void solve(System next){
    a = PVector.sub(R, PVector.mult(next.a,Z));
    b = PVector.sub(next.getPos().mult(2.0), next.a);
  }
  
  System midPoint(System next){
    PVector mid = new PVector();
    mid.add(getPos().mult(0.125));
    mid.add(a.copy().mult(0.375));
    mid.add(b.copy().mult(0.375));
    mid.add(next.getPos().mult(0.125));
    return new System(mid.x, mid.y);
  }
  
  void show(System next){
    float t;
    PVector curve = new PVector();
    for(int i=0; i<=n; i++){
      t = (float)i/n;
      curve.add(getPos().mult((1-t)*(1-t)*(1-t)));
      curve.add(a.copy().mult(3*t*(1-t)*(1-t)));
      curve.add(b.copy().mult(3*t*t*(1-t)));
      curve.add(next.getPos().mult(t*t*t));
      vertex(curve.x, curve.y);
      curve.set(0,0);
    }
  }
  
}

/************************************************************/
Fluid sim;

void setup(){
  size(1280, 700);
  sim = new Fluid(4500);
  noStroke();
}

void draw(){
  float startFrame = millis();
  sim.update();
  float endPhysics = millis();
  
  background(50);
  sim.display();
  
  float endFrame = millis();
  
  String Report = "Display : "+(endFrame-endPhysics)+"ms,"
    +" Physics: "+(endPhysics-startFrame)+"ms,"
    +" FPS: "+(round(frameRate)) + "\n";
  surface.setTitle("GS20119 SPHFLUID - "+Report);
}

void mouseClicked(){
  if(mouseButton==RIGHT) sim.addDrop(mouseX, mouseY);
}
