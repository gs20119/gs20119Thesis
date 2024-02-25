import java.util.LinkedList;
import interfascia.*;

int N;
final float dx=5;
final float dt=0.1;
final float RELAX=1.0;
final int SIZE=500;
final float eps=0.0001;
final float viscosity=2;
final float density=0.1;
float clickX, clickY;
boolean clicked=false;
PVector[] adj = new PVector[8];
Simulation sim;

class Simulation{
   Cell[] Grid;
   ArrayList<Polygon> Bubbles;
   PVector[] boundaryVelocity;
   PVector[] gridVelocity;
   PVector[] projVelocity;
   PVector[] w, b;
   int[] boundaryNormal; 
   boolean[] isBoundary; 
   float[] p, dw;
   PVector[] dp;
     
   Simulation(){
      Grid = new Cell[N*N]; 
      boundaryVelocity = new PVector[N*N];
      boundaryNormal = new int[N*N];
      isBoundary = new boolean[N*N];
      gridVelocity = new PVector[N*N];
      projVelocity = new PVector[N*N];
      w = new PVector[N*N];
      b = new PVector[N*N];
      p = new float[N*N];
      dw = new float[N*N];
      dp = new PVector[N*N];

      Bubbles = new ArrayList<>();
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         Grid[N*i+j] = new Cell(i,j);
         boundaryVelocity[N*i+j] = new PVector();                       
         gridVelocity[N*i+j] = new PVector();
         dp[N*i+j] = new PVector();
      }
      
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         isBoundary[N*i+j]=false;
         if(i==0 || i==N-1 || j==0 || j==N-1){
            isBoundary[N*i+j]=true;
            boundaryVelocity[N*i+j].set(0,0);
            if(i==0) boundaryNormal[N*i+j] = 0;
            else if(j==0) boundaryNormal[N*i+j] = 6;
            else if(i==N-1) boundaryNormal[N*i+j] = 4;
            else boundaryNormal[N*i+j] = 2;
         }
      }
      
   }
   
   void show(){
      loadPixels();
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]) Grid[N*i+j].setColor(200,0,200);
         else{
            PVector v = Grid[N*i+j].getVelocity();
            float p = Grid[N*i+j].getPressure();
            Grid[N*i+j].setColor(round(v.x*10), round(v.y*10), round(-v.x*10));
         }
         Grid[N*i+j].show();
      }
      updatePixels();
      for(Polygon P:Bubbles){
        P.calculate();
        P.show();
      }
   }

   void addBubble(float x, float y){
      Polygon bubble = new Polygon(x,y,20);
      Bubbles.add(bubble);
   }
   
   void solveAdvectionDiffusion(){  
      for(int i=0; i<N; i++) for(int j=0; j<N; j++)
         gridVelocity[N*i+j] = Grid[N*i+j].getVelocity();
      int iter=0; float Error;  
      PVector Sum = new PVector();
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]) continue;
         float I = min(max(i*dx-gridVelocity[N*i+j].x*dt,0),(N-1)*dx);
         float J = min(max(j*dx-gridVelocity[N*i+j].y*dt,0),(N-1)*dx);
         b[N*i+j] = gridVelocity[N*(round(I/dx))+round(J/dx)].copy();
      }
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]){ w[N*i+j] = boundaryVelocity[N*i+j].copy(); continue; }
         if(isBoundary[N*(i+1)+j]) b[N*i+j].add(boundaryVelocity[N*(i+1)+j].copy().mult(viscosity*dt/(dx*dx)));
         if(isBoundary[N*(i-1)+j]) b[N*i+j].add(boundaryVelocity[N*(i-1)+j].copy().mult(viscosity*dt/(dx*dx)));
         if(isBoundary[N*i+(j+1)]) b[N*i+j].add(boundaryVelocity[N*i+(j+1)].copy().mult(viscosity*dt/(dx*dx)));
         if(isBoundary[N*i+(j-1)]) b[N*i+j].add(boundaryVelocity[N*i+(j-1)].copy().mult(viscosity*dt/(dx*dx)));
         w[N*i+j] = b[N*i+j].copy();
      }
      while(true){
         Error=0.0; iter++;
         for(int i=0; i<N; i++) for(int j=0; j<N; j++){
            if(isBoundary[N*i+j]) continue;
            Sum.set(0,0);
            if(!isBoundary[N*(i+1)+j]) Sum.add(w[N*(i+1)+j]);
            if(!isBoundary[N*(i-1)+j]) Sum.add(w[N*(i-1)+j]);
            if(!isBoundary[N*i+(j+1)]) Sum.add(w[N*i+(j+1)]);
            if(!isBoundary[N*i+(j-1)]) Sum.add(w[N*i+(j-1)]);
            Sum.mult(-viscosity*dt/(dx*dx));
            PVector Delta = b[N*i+j].copy().sub(Sum).div(1+4*viscosity*dt/(dx*dx)).sub(w[N*i+j]).mult(RELAX);
            Error += sq(mag(Delta.x, Delta.y));
            w[N*i+j].add(Delta);
         }
         if(sqrt(Error/(N*N))<0.1) break;
      }
   }

   void solvePressure(){
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         p[N*i+j] = Grid[N*i+j].getPressure();
         if(isBoundary[N*i+j]){ dw[N*i+j]=0; continue; }
         dw[N*i+j] = ((w[N*(i+1)+j].x-w[N*(i-1)+j].x)+(w[N*i+(j+1)].y-w[N*i+(j-1)].y))/(2*dx);
         dw[N*i+j] *= density/dt;
      }
      float Error, Sum; int iter=0;
      while(true){
         Error=0.0; iter++;
         for(int i=0; i<N; i++) for(int j=0; j<N; j++){
            Sum=0.0;
            if(isBoundary[N*i+j]){
               int k = boundaryNormal[N*i+j];
               Sum = p[N*(i+floor(adj[k].x))+j+floor(adj[k].y)];
               float Delta = RELAX*(Sum-p[N*i+j]);
               Error += sq(Delta);
               p[N*i+j] += Delta;
            }else{
               Sum += p[N*(i+1)+j]; Sum += p[N*(i-1)+j];
               Sum += p[N*i+(j+1)]; Sum += p[N*i+(j-1)];
               Sum *= 1/(dx*dx);
               float Delta = RELAX*((dw[N*i+j]-Sum)/(-4/(dx*dx))-p[N*i+j]);
               Error += sq(Delta);
               p[N*i+j] += Delta;
            }
         }
         if(sqrt(Error/(N*N))<0.1) break;
      }
   }

   void setValues(){
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]) continue;
         dp[N*i+j].set((p[N*(i+1)+j]-p[N*(i-1)+j])/(2*dx),(p[N*i+(j+1)]-p[N*i+(j-1)])/(2*dx));
         dp[N*i+j].mult(dt/density);
      }
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]) Grid[N*i+j].setVelocity(boundaryVelocity[N*i+j]);
         else Grid[N*i+j].setVelocity(gridVelocity[N*i+j].sub(dp[N*i+j]));
         Grid[N*i+j].setPressure(p[N*i+j]);
      }
   }
   
   void addForces(float x, float y, float X, float Y){
      PVector F = new PVector(X-x, Y-y).normalize().mult(20.0);
      int I = round(x/dx), J = round(y/dx);
      for(int i=0; i<N; i++) for(int j=0; j<N; j++)
        gridVelocity[N*i+j].add(PVector.mult(F,dt*exp(-(sq(I-i)+sq(J-j))/4)));
   }

}

class Cell{
   color dye;
   PVector velocity;
   PVector position;
   float pressure;
   int x, y;
   Cell(int i, int j){
      dye = color(200, 200, 200);
      x = i; y = j;
      position = new PVector(x*dx, y*dx);
      velocity = new PVector(0,0);
      pressure = 1;
   }
   void show(){
      if(pressure<0){ println(x + " " + y + " " + pressure); delay(10); }
      for(int i=round((x)*dx); i<(x+1)*dx; i++)
        for(int j=round((y)*dx); j<(y+1)*dx; j++)
          pixels[SIZE*i+j]=dye;
   }
   void setColor(int x, int y, int z){
      dye = color(x,y,z);
   }
   PVector getVelocity(){
      return velocity.copy();
   }
   void setVelocity(PVector v){
      velocity = v.copy();
   }
   float getPressure(){
      return pressure;
   }
   void setPressure(float p){
      pressure=p;
   }
}


class System{
  PVector a, b; // need to solve
  PVector r, R; 
  float Z; 
  Coordinate P;
  static final int n = 10;
  
  System(float X, float Y){
    a = new PVector();
    b = new PVector();
    r = new PVector();
    R = new PVector();
    P = new Coordinate(X,Y);
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
      for(int i=0; i<vtx; i++){
         if(i==0) r = 200;
         else if(i==vtx-1 || i==1) r = 100;
         else if(i==vtx-2 || i==2) r = 80;
         else if(i==vtx-3 || i==3) r = 75;
         else r = 80;
         //r += random(20)-10;
         sys[i] = new System(x+r*cos(2*PI*i/n+1),y+r*sin(2*PI*i/n+1));
      }
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
      stroke(255,0,0);
      strokeWeight(5);
      for(int i=0; i<vtx; i++){
        PVector pos = sys[i].P.getPos();
        point(pos.x, pos.y);
        stroke(0);
      }
   }
}


class Coordinate{
   float x, y;
   PVector velocity;
   PVector acceleration;
   Coordinate(float X, float Y){
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

void setup(){
   size(500,500);
   N = floor(SIZE/(dx));
   adj[0] = new PVector(1,0);
   adj[1] = new PVector(1,-1);
   adj[2] = new PVector(0,-1);
   adj[3] = new PVector(-1,-1);
   adj[4] = new PVector(-1,0);
   adj[5] = new PVector(-1,1);
   adj[6] = new PVector(0,1);
   adj[7] = new PVector(1,1);
   sim = new Simulation();
   sim.addBubble(250,250);
}

void draw(){
   background(0);
   sim.solveAdvectionDiffusion();
   sim.solvePressure();
   if(!clicked && mousePressed){
      clicked=true;
      clickX=mouseX; clickY=mouseY;
   }
   if(clicked && !mousePressed){
      clicked=false;
      sim.addForces(clickX,clickY,mouseX,mouseY);
   }
   sim.setValues();
   sim.show();
}
