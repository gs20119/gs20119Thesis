import java.util.ArrayList;

class Simulation{
   ArrayList<Integer>[] regions;
   ArrayList<Particle> particles;
   int dx=6, N=4000;
   int I, J;
   
   Simulation(){
     I = width/dx+1; J = height/dx+1; // coord(100,100) <-> size(600,600)
     regions = new ArrayList[I*J];
     for(int i=0; i<I*J; i++) regions[i]=new ArrayList<>();
     particles = new ArrayList<>();
     for(int i=0; i<N; i++){
       float x = random(100);
       float y = random(100);
       particles.add(new Particle(x,y)); 
       particles.get(i).setColor(color(100+x, 100+y, 100-x));
     }
     for(int n=0; n<N; n++){
         PVector r = particles.get(n).getR(); 
         regions[round(r.x)*J+round(r.y)].add(n);
      }
   }
   
   void addParticle(float x, float y, float vx, float vy){
      particles.add(new Particle(x,y,vx,vy)); 
      PVector r = particles.get(N).getR(); 
      particles.get(N).setColor(color(100+x, 100+y, 100-x));
      regions[round(r.x)*J+round(r.y)].add(N); N++;
   }
   
   float W(PVector r){
      float c = 0.266;
      float R = r.mag();
      float h = 50; // smoothing distance
      float q = R/h;
      if(q<=0.5) return c*6*(pow(q,3)-pow(q,2))+1;
      else if(q<=1) return c*2*(pow(1-q,3));
      else return 0;
   }
   
   PVector dW(PVector r){
      float c = 0.266;
      float R = r.mag();
      float h = 50;
      float q = R/h;
      if(q<=0.5) return r.mult(c*6*(3*R-2*h)/pow(h,3));
      else if(q<=1) return r.mult(-c*6*pow(h-R,2)/(pow(h,3)*R));
      else return new PVector();
   }
   
   void Simulate(){
      solveNonpForce();
      move();
      solveDensityPressure();
      solvePForce();
      move();
      show();
   }
   
   void solveNonpForce(){
      float nu = 0.00;
      for(int n=0; n<N; n++){
         PVector nr = particles.get(n).getR();
         PVector nv = particles.get(n).getV();
         float nm = particles.get(n).getM();
         float np = particles.get(n).getP();
         float nrho = particles.get(n).getRho();
         PVector acc = new PVector(0,0);
         int ii = round(nr.x), jj = round(nr.y); // (ii,jj) : center
         for(int i=max(ii-5,0); i<=min(ii+5,I-1); i++)
           for(int j=max(jj-5,0); j<=min(jj+5,J-1); j++){ // (i,j) : neigborhood
             if((ii-i)*(ii-i)+(jj-j)*(jj-j)>26) continue;
             for(int m : regions[i*J+j]){ // particles.get(m]
                PVector mr = particles.get(m).getR();
                PVector mv = particles.get(m).getV();
                float mm = particles.get(m).getM();
                float mp = particles.get(m).getP();
                float mrho = particles.get(m).getRho();
                PVector rr = PVector.sub(nr,mr);
                PVector vv = PVector.sub(nv,mv);
                acc.add(vv.mult(-W(rr)*nm*nu*mm/mrho));
                //acc.add(dW(rr).mult(nu*8*mm*(rr.x*vv.x+rr.y*vv.y)/(mrho*(pow(rr.mag(),2)+0.01))));
             }
         }
         acc.add(new PVector(0,0)); // external
         acc.div(nm);
         particles.get(n).setA(acc);
      }
   }
   
   void solveDensityPressure(){
      float k = 1;
      float dt = 0.1;
      float density0 = particles.get(0).getM()*N/(I*J);
      for(int n=0; n<N; n++){
         PVector nr = particles.get(n).getR();
         PVector nv = particles.get(n).getV();
         float density = 0;
         int ii = round(nr.x), jj = round(nr.y); // (ii,jj) : center
         for(int i=max(ii-5,0); i<=min(ii+5,I-1); i++)
           for(int j=max(jj-5,0); j<=min(jj+5,J-1); j++){ // (i,j) : neigborhood
             if((ii-i)*(ii-i)+(jj-j)*(jj-j)>25.1) continue;
             for(int m : regions[i*J+j]){ 
                PVector mr = particles.get(m).getR();
                PVector mv = particles.get(m).getV();
                float mm = particles.get(m).getM();
                PVector rr = PVector.sub(nr,mr);
                PVector vv = PVector.sub(nv,mv);
                density += mm*W(rr);
                density += dt*mm*(vv.x*dW(rr).x+vv.y*dW(rr).y);
             }
         }
         //println(density);
         println(k*(density-density0));
         particles.get(n).setRho(density);
         particles.get(n).setP(k*(density/density0-1));
      }
   }
   
   
   
   void solvePForce(){
      float dt = 0.1;
      for(int n=0; n<N; n++){
         PVector nr = particles.get(n).getR();
         PVector nv = particles.get(n).getV();
         float nm = particles.get(n).getM();
         float np = particles.get(n).getP();
         float nrho = particles.get(n).getRho();
         PVector acc = new PVector();
         int ii = round(nr.x), jj = round(nr.y); // (ii,jj) : center
         for(int i=max(ii-5,0); i<=min(ii+5,I-1); i++)
           for(int j=max(jj-5,0); j<=min(jj+5,J-1); j++){ // (i,j) : neigborhood
             if(ii==i && jj==j) continue;
             for(int m : regions[i*J+j]){ // particles.get(m]
                PVector mr = particles.get(m).getR();
                float mm = particles.get(m).getM();
                float mp = particles.get(m).getP();
                float mrho = particles.get(m).getRho();
                PVector rr = PVector.sub(nr,mr);
                acc.add(dW(rr).mult(-mm*(mp/(mrho*mrho)+np/(nrho*nrho))));
             }
         }
         acc.div(nm);
         particles.get(n).setA(acc);
      }
   }
   
   void move(){
      float dt = 0.1;
      for(int n=0; n<N; n++){
         PVector nr = particles.get(n).getR();
         PVector nv = particles.get(n).getV();
         PVector na = particles.get(n).getA();
         int prevReg = round(nr.x)*J+round(nr.y);
         nv.add(PVector.mult(na,dt/2.0));
         nr.add(PVector.mult(nv,dt));
         nv.add(PVector.mult(na,dt/2.0));
         if(nr.x<0 || nr.x>width/dx){ nr.x = constrain(nr.x, 0.001, width/dx-0.001); nv.x *= -0.1; }
         if(nr.y<0 || nr.y>height/dx){ nr.y = constrain(nr.y, 0.001, height/dx-0.001); nv.y *= -0.1; }
         int nextReg = round(nr.x)*J+round(nr.y);
         if(prevReg != nextReg){
            regions[prevReg].remove(regions[prevReg].indexOf(n));
            regions[nextReg].add(n);
         }
         particles.get(n).setR(nr);
         particles.get(n).setV(nv);
         particles.get(n).setA(new PVector());
      }
   }
   
   void show(){
      for(int n=0; n<N; n++){
         PVector nr = particles.get(n).getR();
         noStroke();
         fill(particles.get(n).getColor());
         circle(nr.x*dx, nr.y*dx, 10);
      }
   }
}

class Particle{
   float p=0;
   float rho=2000*0.3/10000;
   float m;
   PVector v, r, a; // r : in real size
   PVector pressGrad;
   color c;
   Particle(float x, float y, float vx, float vy){
      m = 0.3;
      v = new PVector(vx,vy);
      r = new PVector(x,y);
      a = new PVector();
      c = color(200);
   }
   Particle(float x, float y){
      m = 0.3;
      v = new PVector();
      r = new PVector(x,y);
      a = new PVector();
      c = color(200);
   }
   void setRho(float density){ rho = density; }
   void setP(float pressure){ p = pressure; }
   void setA(PVector acc){ a = acc.copy(); }
   void setV(PVector velocity){ v = velocity.copy(); }
   void setR(PVector position){ r = position.copy(); }
   void setColor(color dye){ c = dye; }
   
   PVector getR(){ return r.copy(); }
   PVector getV(){ return v.copy(); }
   PVector getA(){ return a.copy(); }
   float getRho(){ return rho; }
   float getP(){ return p; }
   float getM(){ return m; }
   color getColor(){ return c; }
   
}

Simulation sim;
int iter=0;
color dye = color(100,200,100);

void setup(){
  size(600,600);
  sim = new Simulation();
}

void draw(){
  println(iter++);
  background(0);
  sim.Simulate();
}
