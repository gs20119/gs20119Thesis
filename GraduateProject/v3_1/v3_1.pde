import java.util.LinkedList;
import interfascia.*;

int N;
final float dx=5, dt=1;
final float RELAX=1;
final int SIZE=500;
final float eps=0.0001;
final float viscosity=1;
PVector[] adj = new PVector[8];
Simulation sim;

class Simulation{
   Cell[] Grid;
   ArrayList<Polygon> Bubbles;
   PVector[] boundaryVelocity;
   PVector[] gridVelocity;
   PVector[] projVelocity;
   PVector[] x, b;
   int[] boundaryNormal; 
   boolean[] isBoundary; 
   float[][] A;
     
   Simulation(){
      Grid = new Cell[N*N]; 
      boundaryVelocity = new PVector[N*N];
      boundaryNormal = new int[N*N];
      isBoundary = new boolean[N*N];
      gridVelocity = new PVector[N*N];
      projVelocity = new PVector[N*N];
      x = new PVector[N*N];
      b = new PVector[N*N];
      A = new float[N*N][N*N];

      Bubbles = new ArrayList<>();
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         Grid[N*i+j] = new Cell(i,j);
         boundaryVelocity[N*i+j] = new PVector();
         gridVelocity[N*i+j] = new PVector();
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
         else Grid[N*i+j].setColor(200,200,200);
         Grid[N*i+j].show();
      }
      for(Polygon P:Bubbles) P.show();
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
         float I = min(max(i*dx-gridVelocity[N*i+j].x*dt,0),(N-1)*dx);
         float J = min(max(j*dx-gridVelocity[N*i+j].y*dt,0),(N-1)*dx);
         b[N*i+j] = gridVelocity[N*(round(I/dx))+round(J/dx)].copy();
      }
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         if(isBoundary[N*i+j]) continue;
         if(isBoundary[N*(i+1)+j]) b[N*i+j].add(boundaryVelocity[N*(i+1)+j].copy().mult(viscosity*dt/(dx*dx)));
         if(isBoundary[N*(i-1)+j]) b[N*i+j].add(boundaryVelocity[N*(i-1)+j].copy().mult(viscosity*dt/(dx*dx)));
         if(isBoundary[N*i+(j+1)]) b[N*i+j].add(boundaryVelocity[N*i+(j+1)].copy().mult(viscosity*dt/(dx*dx)));
         if(isBoundary[N*i+(j-1)]) b[N*i+j].add(boundaryVelocity[N*i+(j-1)].copy().mult(viscosity*dt/(dx*dx)));
         x[N*i+j] = b[N*i+j].copy();
      }
      while(true){
         Error=0.0; iter++;
         for(int i=0; i<N; i++) for(int j=0; j<N; j++){
            if(isBoundary[N*i+j]) continue;
            Sum.set(0,0);
            if(!isBoundary[N*(i+1)+j]) Sum.add(x[N*(i+1)+j]); 
            if(!isBoundary[N*(i-1)+j]) Sum.add(x[N*(i-1)+j]);
            if(!isBoundary[N*i+(j+1)]) Sum.add(x[N*i+(j+1)]);
            if(!isBoundary[N*i+(j-1)]) Sum.add(x[N*i+(j-1)]);
            Sum.mult(-viscosity*dt/(dx*dx));
            PVector Delta = b[N*i+j].copy().sub(Sum).div(1+4*viscosity*dt/(dx*dx)).sub(x[N*i+j]).mult(RELAX);
            Error += sq(mag(Delta.x, Delta.y));
            x[N*i+j].add(Delta);
         }
         if(sqrt(Error)<0.01) break;
      }println(iter);
   }

   void solvePressure(){
      
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
      pressure = 1.0;
   }
   void show(){
      stroke(dye.x, dye.y, dye.z);
      strokeWeight(dx*1.25);
      point(x*dx,y*dx);
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
   void setVelocity(float x, float y){
      velocity.set(x,y);
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
   N = floor(SIZE/dx)+1;
   adj[0] = new PVector(1,0);
   adj[1] = new PVector(1,-1);
   adj[2] = new PVector(0,-1);
   adj[3] = new PVector(-1,-1);
   adj[4] = new PVector(-1,0);
   adj[5] = new PVector(-1,1);
   adj[6] = new PVector(0,1);
   adj[7] = new PVector(1,1);
   sim = new Simulation();
   sim.addBubble(300,300);
   sim.addBubble(200,200);
}

void draw(){
   background(0);
   sim.getBoundary();
   sim.solveAdvectionDiffusion();
   sim.show();
}
