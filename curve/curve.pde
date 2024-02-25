import java.util.LinkedList;

final float dx=1;
final float SCALE=5;

class ColorPoint{
   float x, y;
   PVector velocity;
   PVector acceleration;
   ColorPoint(float X, float Y){
      x = X; y = Y;
      velocity = new PVector();
      acceleration = new PVector();
   }
   PVector getPos(){
      return new PVector(x,y);
   }
   void setPos(float X, float Y){
      x = X; y = Y;
   }
   PVector getVel(){
      return velocity.copy();
   }
}

class System{
  PVector a, b; // need to solve
  PVector r, R; 
  float Z; 
  ColorPoint P;
  static final int n = 10;
  
  System(float X, float Y){
    a = new PVector();
    b = new PVector();
    r = new PVector();
    R = new PVector();
    P = new ColorPoint(X,Y);
  }
  
  PVector getPos(){
    return P.getPos();  
  }
  
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

class Polygon{
   int vtx;
   System[] sys;
   PVector dye;
   
   Polygon(float x, float y, int n){
      sys = new System[n];
      vtx = n;
      dye = new PVector(100,150,200);
      float r = 100;
      for(int i=1; i<vtx-2; i++){
         r = 100;
         r += random(40)-20;
         sys[i] = new System(x+r*cos(2*PI*i/n+1),y+r*sin(2*PI*i/n+1));
      }
      float r1 = 100+random(60)-30;
      float r2 = 100+random(60)-30;
      float r3 = (r1+r2)*0.5;
      sys[vtx-2] = new System(x+r*cos(2*PI*(vtx-2)/n+1),y+r*sin(2*PI*(vtx-2)/n+1));
      sys[vtx-1] = new System(x+r*cos(2*PI*(vtx-1)/n+1),y+r*sin(2*PI*(vtx-1)/n+1));
      sys[0] = new System(x+r*cos(2*PI*0/n+1),y+r*sin(2*PI*0/n+1));
   }
   void calculate(){
      sys[0].setValue(sys[1]);
      for(int i=1; i<vtx-1; i++)
        sys[i].setValue(sys[i-1], sys[i+1]);
      sys[vtx-1].setValue(sys[vtx-2], sys[0], 0);
      sys[vtx-1].solve(sys[0], 0);
      for(int i=vtx-2; i>=0; i--)
        sys[i].solve(sys[i+1]);
   }
   void show(){
      noStroke();
      smooth();
      fill(dye.x, dye.y, dye.z);
      beginShape();
      for(int i=0; i<vtx-1; i++)
        sys[i].show(sys[i+1]);
      sys[vtx-1].show(sys[0]);
      endShape();
      beginShape();
      for(int i=0; i<vtx; i++){
        PVector pos = sys[i].P.getPos();
        vertex(pos.x+300, pos.y);
      }
      vertex(sys[0].P.getPos().x+300, sys[0].P.getPos().y);
      endShape();
   }
}

Polygon p;

void setup(){
  size(1000,500);
  background(255);
  p = new Polygon(250, 250, 15);
}

void draw(){
  p.calculate();
  p.show();
}
