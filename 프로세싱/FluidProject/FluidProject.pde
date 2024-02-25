
class Particle{
  PVector position;
  PVector velocity;
  PVector acceleration;
  float lifespan;
  float mass;
  float density;
  float pressure;
  float t=0, dt=0.01;
  float scale;
  
  Particle(PVector p, float sc){
     mass = 0.05;
     acceleration = new PVector(0, 0);
     velocity = new PVector(0, 0);
     position = p.copy();
     lifespan = 255.0;
     density = 0;
     pressure = 0;
     scale = sc;
  }
  
  void run(){
    update();
    display();
    t += dt;
  }
  
  void update(){
    velocity.add(acceleration.copy().mult(0.5*dt));
    position.add(velocity.copy().mult(dt));
    if(position.x<0){ position.x=0; velocity.x *= -1; }
    if(position.x>scale){ position.x=scale; velocity.x *= -1; }
    if(position.y<0){ position.y=0; velocity.y *= -1; }
    if(position.y>scale){ position.y=scale; velocity.y *= -1; }
    velocity.add(acceleration.copy().mult(0.5*dt));
  }
  
  void display(){
    stroke(255, 25*density);
    fill(255, 25*density);
    ellipse(position.x*width/scale, position.y*height/scale, 20, 20);
  }
  
}

class ParticleSystem{
   ArrayList<Particle> particles;
   PVector[][] pipj;
   float h=0.15, k=0.2, n=1, nu=0.7;
   float scale=10;
   
   ParticleSystem(int count){
      pipj = new PVector[count][count];
      particles = new ArrayList<Particle>();
      for(int i=0; i<count; i++)
      particles.add(new Particle(new PVector(random(0,scale),random(scale*1/3,scale)), scale));
   }
   
   void run(){
      for(int i=0; i<particles.size(); i++) particles.get(i).run();
      for(int i=0; i<particles.size(); i++){
         for(int j=0; j<particles.size(); j++){
           pipj[i][j] = particles.get(i).position.copy();
           pipj[i][j].sub(particles.get(j).position);
         }
      }getDP(); getAcc();
      System.out.printf("%.2f\n", particles.get(0).t);
   }
   
   float W(PVector p1p2){
     float dist = p1p2.mag();
     float w = exp(-pow(dist/h,2))/pow(h*sqrt(PI),3);
     return w;
  }
  
  PVector gradW(PVector p1p2){
     float dist = p1p2.mag();
     float n = -2*exp(-pow(dist/h,2))/(pow(h,5)*pow(PI,1.5));
     PVector gw = p1p2.copy().mult(n);
     return gw;
  }
   
   void getDP(){
      for(int i=0; i<particles.size(); i++){
         float density=0;
         Particle pi = particles.get(i);
         for(int j=0; j<particles.size(); j++){
            Particle pj = particles.get(j);
            density += pj.mass * W(pipj[i][j]);
         }
         pi.density = density;
         pi.pressure = k*pow(density,1+1/n);
      }
   }
   
   void getAcc(){
     for(int i=0; i<particles.size(); i++){
       PVector force = new PVector(0,0);
       Particle pi = particles.get(i);
       for(int j=0; j<particles.size(); j++){
          if(i==j) continue;
          Particle pj = particles.get(j);
          float pressure = pi.pressure/pow(pi.density,2) + pj.pressure/pow(pj.density,2);
          force.add(gradW(pipj[i][j]).copy().mult(-pj.mass*pressure));
       }
       force.add(pi.velocity.copy().mult(-nu));
       if(pi.position.x>scale*0.4 && pi.position.x<scale*0.6)
         force.add(new PVector(0,-2*pi.position.y/scale));
       else if(pi.position.y>scale*0.8) force.add(new PVector(-3*(pi.position.x-scale*0.5)/scale,0));
       force.add(new PVector(0,0.5));
       pi.acceleration = force;
     }
   }
}

ParticleSystem ps;

void setup(){
  size(640, 640);
  ps = new ParticleSystem(1500);
  frameRate(120);
}

void draw(){
  background(0);
  ps.run();
}
