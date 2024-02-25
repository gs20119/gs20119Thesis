
int N;
float dx=5, dt=1;
final int SIZE=500;
Cell[][] Grid;

class Cell{
   PVector dye;
   PVector velocity;
   float pressure;
   PVector position;
   Cell(int i, int j){
      dye = new PVector(200, 200, 200);
      position = new PVector(i,j);
      velocity = new PVector(0,0);
      pressure = 1.0;
   }
   void show(){
      stroke(dye.x, dye.y, dye.z);
      strokeWeight(dx*1.25);
      point(position.x*dx,position.y*dx);
   }
}

class Polygon{
   int num_vertex;
   ArrayList<Integer> points;
   ArrayList<Coordinate> coords;
   Polygon(float x, float y, float r, int n){
      points = new ArrayList<>();
      coords = new ArrayList<>();
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
         PVector cd = coords.get(points.get(i)).position;
         point(cd.x, cd.y);
      }
   }
}

class Coordinate{
   PVector position;
   PVector velocity;
   PVector acceleration;
   Coordinate(float x, float y){
      position = new PVector(x,y);
      velocity = new PVector();
      acceleration = new PVector();
   }
}

void setup(){
   size(500,500);
   N = floor(SIZE/dx)+1;
   Grid = new Cell[N][N]; 
   for(int i=0; i<N; i++) for(int j=0; j<N; j++) 
      Grid[i][j] = new Cell(i,j);
}

void draw(){
   background(0);
   for(int i=0; i<N; i++) for(int j=0; j<N; j++)
      Grid[i][j].show();
   Polygon ten = new Polygon(300, 300, 100, 10);
   ten.show();
}
