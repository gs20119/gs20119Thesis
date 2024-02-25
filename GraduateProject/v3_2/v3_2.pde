import java.util.LinkedList;
import interfascia.*;

int N;
final float dx=1, dt=0.1;
final float RELAX=0.9;
final float SCALE=5;
final int SIZE=500;
final float eps=0.0001;
final float viscosity=2;
final float density=0.1;
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
   }

   void getBoundary(){
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
      for(Polygon P:Bubbles)
         for(int i=0; i<N; i++) for(int j=0; j<N; j++){
            int edge = Grid[N*i+j].isBoundary(P);
            if(edge!=-1){
               isBoundary[N*i+j]=true;
               boundaryVelocity[N*i+j] = Grid[N*i+j].getBoundaryVelocity(P,edge);
               boundaryNormal[N*i+j] = Grid[N*i+j].getBoundaryNormal(P,edge);
            }
         }
   }

   void show(){
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]) Grid[N*i+j].setColor(200,0,200);
         else{
            PVector v = Grid[N*i+j].getVelocity();
            float p = Grid[N*i+j].getPressure();
            Grid[N*i+j].setColor(round(v.x*10), round(-v.x*10), 0);//round(p));
         }
         Grid[N*i+j].show();
      }
      for(Polygon P:Bubbles) P.show();
      stroke(0); strokeWeight(3);
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]){
            int k = boundaryNormal[N*i+j];
            //line(i*dx, j*dx, (i+floor(adj[k].x))*dx, (j+floor(adj[k].y))*dx);
         }
      }
   }

   void addBubble(float x, float y){
      Polygon bubble = new Polygon(x,y,65,10);
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
         else Grid[N*i+j].setVelocity(gridVelocity[N*i+j].copy().sub(dp[N*i+j]));
         Grid[N*i+j].setPressure(p[N*i+j]);
      }
   }

}

class Cell{
   PVector dye;
   PVector velocity;
   PVector position;
   float pressure;
   int x, y;
   Cell(int i, int j){
      dye = new PVector(200, 200, 200);
      x = i; y = j;
      position = new PVector(x*dx, y*dx);
      velocity = new PVector(0,0);
      pressure = 1;
   }
   void show(){
      if(pressure<0){ println(x + " " + y + " " + pressure); delay(10); }
      stroke(dye.x, dye.y, dye.z);
      strokeWeight(dx*1.25*SCALE);
      point(x*dx*SCALE,y*dx*SCALE);
   }
   int isBoundary(Polygon P){
      if(x==0 || x==N-1 || y==0 || y==N-1) return -1;
      for(int i=0; i<8; i++){
         int edge = P.isIntersect(position, position.copy().add(adj[i].copy().mult(dx)));
         if(edge!=-1) return edge;
      }return -1;
   }
   PVector getBoundaryVelocity(Polygon P, int edge){
      return P.edgeVelocity(position, edge);
   }
   int getBoundaryNormal(Polygon P, int edge){
      return P.edgeNormal(position, edge);
   }
   void setColor(int x, int y, int z){
      dye.set(x,y,z);
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

int CCW(PVector a, PVector b, PVector c){
   float op = a.x*b.y+b.x*c.y+c.x*a.y;
   op -= a.y*b.x+b.y*c.x+c.y*a.x;
   if(op>0) return 1;
   else if(op==0) return 0;
   else return -1;
}

class Polygon{
   int num_vertex;
   LinkedList<Integer> points;
   LinkedList<Coordinate> coords;
   Polygon(float x, float y, float r, int n){
      points = new LinkedList<>();
      coords = new LinkedList<>();
      num_vertex = n;
      for(int i=0; i<n; i++){
         points.add(i);
         coords.add(new Coordinate(x+r*cos(2*PI*i/n+1),y+r*sin(2*PI*i/n+1)));
      }
   }
   void show(){
      stroke(0);
      strokeWeight(dx*1.25);
      for(int i=0; i<num_vertex; i++){
         PVector cd = coords.get(points.get(i)).getPos();
         point(cd.x, cd.y);
      }
   }
   int isIntersect(PVector a, PVector b){
      for(int i=0; i<num_vertex; i++){
         PVector c = coords.get(points.get(i)).getPos();
         PVector d = coords.get(points.get((i+1)%num_vertex)).getPos();
         int ab = CCW(a,b,c)*CCW(a,b,d);
         int cd = CCW(c,d,a)*CCW(c,d,b);
         if(ab<=0 && cd<=0) return i;
      }return -1;
   }
   PVector edgeVelocity(PVector pos, int edge){
      PVector v1 = coords.get(points.get(edge)).getPos();
      PVector v2 = coords.get(points.get((edge+1)%num_vertex)).getPos();
      float a = v2.y-v1.y;
      float b = v1.x-v2.x;
      float c = v1.x*v2.y-v2.x*v1.y;
      float d = (v1.x-v2.x)*pos.x+(v1.y-v2.y)*pos.y;
      float x = (a*c+b*d)/(a*a+b*b);
      PVector V1 = coords.get(points.get(edge)).getVel();
      PVector V2 = coords.get(points.get((edge+1)%num_vertex)).getVel();
      float m = (v1.x-x)/(v1.x-v2.x+eps);
      float n = (x-v2.x)/(v1.x-v2.x+eps);
      return new PVector((V1.x*n+V2.x*m),(V1.y*n+V2.y*m));
   }
   int edgeNormal(PVector pos, int edge){
      PVector v1 = coords.get(points.get(edge)).getPos();
      PVector v2 = coords.get(points.get((edge+1)%num_vertex)).getPos();
      float slope = (v1.x-v2.x)/(v2.y-v1.y+eps);
      int k;
      if(-0.5<slope && slope<=0.5) k=0;
      else if(-2.0<slope && slope<=-0.5) k=1;
      else if(slope<=-2.0 || slope>2.0) k=2;
      else k=3;
      PVector npos = pos.copy().add(adj[k].copy().mult(dx));
      int ab = CCW(pos,npos,v1)*CCW(pos,npos,v2);
      int cd = CCW(v1,v2,pos)*CCW(v1,v2,npos);
      if(ab<=0 && cd<=0) k+=4;
      return k;
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
   N = floor(SIZE/(dx*SCALE))+1;
   adj[0] = new PVector(1,0);
   adj[1] = new PVector(1,-1);
   adj[2] = new PVector(0,-1);
   adj[3] = new PVector(-1,-1);
   adj[4] = new PVector(-1,0);
   adj[5] = new PVector(-1,1);
   adj[6] = new PVector(0,1);
   adj[7] = new PVector(1,1);
   sim = new Simulation();
}

void draw(){
   background(0);
   sim.getBoundary();
   sim.solveAdvectionDiffusion();
   sim.solvePressure();
   sim.setValues();
   sim.show();
}
