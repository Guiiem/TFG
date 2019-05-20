void setup() {
  // put your setup code here, to run once:
Serial.begin(9600);


}

void loop() {
  // put your main code here, to run repeatedly:
int LecturaV = analogRead(A0);

float voltage = LecturaV * (5.0/1023.0);

Serial.print("Voltatge:");
Serial.println(voltage);
Serial.print("Analogic:");
Serial.println(LecturaV);

delay(500);

}
