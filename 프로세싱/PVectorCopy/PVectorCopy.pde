
void setup(){
  PVector a = new PVector(10,10);
  PVector b = new PVector(4,1);
  PVector c = a;
  c.sub(b); 
  System.out.println(c);
  System.out.println(a);
  
  PVector d = new PVector(10,10);
  PVector e = new PVector(4,1);
  PVector f = d.copy();
  f.sub(e);
  System.out.println(f);
  System.out.println(d);
}

void draw(){
  
}
