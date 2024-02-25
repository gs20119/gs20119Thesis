import java.util.LinkedList;
import interfascia.*;

int N;
float dx=5, dt=1;
final int SIZE=500;
PVector[] adj = new PVector[4];
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
   Cell[][] Grid;
   ArrayList<Polygon> Bubbles;
   Vector 
   Simulation(){
      Grid = new Cell[N][N]; 
      Bubbles = new ArrayList<>();
      for(int i=0; i<N; i++) for(int j=0; j<N; j++) 
         Grid[i][j] = new Cell(i,j);
      
   }
   void getBoundary(){
      for(Polygon P:Bubbles)
         for(int i=0; i<N; i++) for(int j=0; j<N; j++){
            int temp = Grid[i][j].isBoundary(P);
            if(temp!=-1) Grid[i][j].setColor(200,0,200);
         }
   }
   void show(){
      for(int i=0; i<N; i++) for(int j=0; j<N; j++)
      Grid[i][j].show();
      for(Polygon P:Bubbles) P.show();
   }

   void addBubble(float x, float y){
      Polygon bubble = new Polygon(x,y,50,20);
      Bubbles.add(bubble);
   }
   
   void solveAdvectionDiffusion(){
      
   }
}

class Cell{
   PVector dye;
   PVector velocity;
   float pressure;
   int x, y;
   Cell(int i, int j){
      dye = new PVector(200, 200, 200);
      x = i; y = j;
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
      PVector org = new PVector(x*dx, y*dx);
      for(int i=0; i<4; i++){
         int temp = P.isIntersect(org,org.copy().add(adj[i]));
         if(temp!=-1) return temp;
      }return -1;
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
}


void setup(){
   size(500,500);
   N = floor(SIZE/dx)+1;
   adj[0] = new PVector(dx,0);
   adj[1] = new PVector(-dx,0);
   adj[2] = new PVector(0,dx);
   adj[3] = new PVector(0,-dx);
   sim = new Simulation();
   sim.addBubble(300,300);
   sim.addBubble(200,200);
}

void draw(){
   background(0);
   sim.getBoundary();
   sim.show();
}
