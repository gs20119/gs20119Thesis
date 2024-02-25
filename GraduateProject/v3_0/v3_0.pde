import java.util.LinkedList;
import interfascia.*;

int N;
final float dx=5, dt=1;
final int SIZE=500;
final float eps=0.0001;
PVector[] adj = new PVector[8];
Simulation sim;

class Vector{
   float[] val;
   Vector(int n){
      val = new float[n];
   }
}

class Matrix{
   float[][] val;
   Matrix(int m, int n){
      val = new float[m][n]; 
   }
}

class Simulation{
   Cell[] Grid;
   ArrayList<Polygon> Bubbles;
   PVector[] boundaryVelocity;
   int[] boundaryNormal; 
   boolean[] isBoundary; 
     
   Simulation(){
      Grid = new Cell[N*N]; 
      boundaryVelocity = new PVector[N*N];
      boundaryNormal = new int[N*N];
      isBoundary = new boolean[N*N];
      Bubbles = new ArrayList<>();
      for(int i=0; i<N; i++) for(int j=0; j<N; j++){
         Grid[N*i+j] = new Cell(i,j);
         boundaryVelocity[N*i+j] = new PVector();
      }
   }
   void getBoundary(){
      for(Polygon P:Bubbles)
         for(int i=0; i<N; i++) for(int j=0; j<N; j++){
            isBoundary[N*i+j]=false;
            int edge = Grid[N*i+j].isBoundary(P);
            if(edge!=-1){
               isBoundary[N*i+j]=true;
               boundaryVelocity[N*i+j] = Grid[N*i+j].getBoundaryVelocity(P,edge);
               boundaryNormal[N*i+j] = Grid[N*i+j].getBoundaryNormal(P,edge);
               Grid[N*i+j].setColor(200,0,200);
            }
         }
   }
   void show(){
      for(int i=0; i<N; i++) for(int j=0; j<N; j++)
      Grid[N*i+j].show();
      for(Polygon P:Bubbles) P.show();
   }

   void addBubble(float x, float y){
      Polygon bubble = new Polygon(x,y,50,20);
      Bubbles.add(bubble);
   }
   
   void solveAdvectionDiffusion(){

      // get Matrix A and Vector b
      // solve Aw=b
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
         int edge = P.isIntersect(position, position.copy().add(adj[i]));
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
      PVector npos = pos.copy().add(adj[k]);
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
   adj[0] = new PVector(dx,0);
   adj[1] = new PVector(dx,-dx);
   adj[2] = new PVector(0,-dx);
   adj[3] = new PVector(-dx,-dx);
   adj[4] = new PVector(-dx,0);
   adj[5] = new PVector(-dx,dx);
   adj[6] = new PVector(dx,0);
   adj[7] = new PVector(dx,dx);
   sim = new Simulation();
   sim.addBubble(300,300);
   sim.addBubble(200,200);
}

void draw(){
   background(0);
   sim.getBoundary();
   sim.show();
}
