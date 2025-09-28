#include <Arduino.h>
#include <Wire.h>
#include <LiquidCrystal_I2C.h>
#include <SD.h>
#include <math.h>
#include <SoftwareSerial.h>

SoftwareSerial Mon(2, 4);

// Helper functions to print to both MATLAB Serial and HW-598 Serial
void printBoth(const String& s){ Serial.print(s); Mon.print(s); }
void printlnBoth(const String& s){ Serial.println(s); Mon.println(s); }
void printBothFloat(float v, int p){ Serial.print(v, p); Mon.print(v, p); }


// ---------------- LCD ----------------
LiquidCrystal_I2C lcd(0x27, 16, 2);

// ---------------- Keypad (TTP229) ----------------
const uint8_t KP_SCL_PIN = 7;
const uint8_t KP_SDO_PIN = 8;

// ---------------- SD ----------------
const int SD_CS_PIN = 10;
File dataFile;
String fileName = "";

// ---------------- Stages ----------------
enum Stage : uint8_t { ASK_D, ASK_d };
Stage stage = ASK_D;

// ---------------- Inputs & state ----------------
uint8_t lastKey = 0;
String currentInput = "";
// values: [0]=D, [1]=d, [2]=P1, [3]=P2
float values[4] = {0};

// ---------------- Sweep control ----------------
bool sweeping = false;
float U0 = 0.0f, U0_start = 0.0f, U0_end = 0.0f, U0_step = 0.1f;

// ---------------- Prompts ----------------
const char* const prompts[2] = {
  "Pipe D (mm):",
  "Orifice d (mm):"
};

// ======================= Keypad Scan =======================
uint8_t scanKeypad() {
  uint8_t key = 0;
  for (uint8_t i = 1; i <= 16; ++i) {
    digitalWrite(KP_SCL_PIN, LOW);
    delayMicroseconds(2);
    if (!digitalRead(KP_SDO_PIN)) key = i;
    digitalWrite(KP_SCL_PIN, HIGH);
    delayMicroseconds(2);
  }
  return key;
}

// ======================= File Helpers =======================
String getNextFileName(bool createNew) {
  uint16_t fileNum = 1;
  String fname;
  while (true) {
    char buf[20];
    snprintf(buf, sizeof(buf), "flow_%03u.csv", fileNum);
    fname = String(buf);
    if (!SD.exists(fname.c_str())) {
      if (createNew) return fname;                 // new: first free number
      else return (fileNum == 1) ? "" : String("flow_") + (fileNum-1 < 100 ? (fileNum-1 < 10 ? "00" : "0") : "") + String(fileNum-1) + ".csv";
    }
    fileNum++;
  }
}

void askFileMode() {
  lcd.clear();
  lcd.setCursor(0,0); lcd.print("1: New file?");
  lcd.setCursor(0,1); lcd.print("2: Old file?");

  while (true) {
    uint8_t k = scanKeypad();
    if (k && k != lastKey) {
      lastKey = k;
      if (k == 1 || k == 2) {
        bool createNew = (k == 1);
        fileName = getNextFileName(createNew);
        if (fileName == "") {
          lcd.clear(); lcd.print("No old files!");
          delay(1500);
          // re-ask
          lcd.clear();
          lcd.setCursor(0,0); lcd.print("1: New file?");
          lcd.setCursor(0,1); lcd.print("2: Old file?");
          continue;
        }
        lcd.clear();
        lcd.print("File: "); lcd.print(fileName);
        delay(1200);

        dataFile = SD.open(fileName.c_str(), FILE_WRITE);
        if (!dataFile) {
          lcd.clear(); lcd.print("File Error!");
          while (1) {}
        }
        // CSV header (only if new/empty)
        if (dataFile.size() == 0) {
          dataFile.println("D_mm,d_mm,U0_mps,P1_Pa,P2_Pa,Qv_m3ps");
          dataFile.flush();
        }
        break;
      }
    }
    delay(40);
  }
}

// ======================= Prompt User =======================
void promptUser() {
  lcd.clear();
  lcd.setCursor(0,0);
  lcd.print(prompts[stage]);   // use global 'stage'
  lcd.setCursor(0,1);
  lcd.print("Type & Enter");
  currentInput = "";
}

// ======================= Flow Calculation =======================
float computeQv(float D_mm, float d_mm, float p1, float p2) {
  const float rho = 997.0f;
  const float epsilon = 1.0f;
  const float C_d = 0.6047f;
  const float PI_4 = PI / 4.0f;

  float delta_p = p1 - p2;
  if (delta_p <= 0) return NAN;

  float D = D_mm / 1000.0f;
  float d = d_mm / 1000.0f;

  float beta = d / D;

  float q_m = C_d / sqrtf(1.0f - powf(beta, 4.0f))
             * epsilon
             * PI_4 * powf(d, 2.0f)
             * sqrtf(2.0f * delta_p * rho);

  float q_v = q_m / rho;

  // Calibration (your earlier linear corrections)
  if (D_mm == 100.0f && d_mm == 50.0f) {
    q_v = 1.027562020238f * q_v + 4.467897137161e-05f;
  } else if (D_mm == 150.0f && d_mm == 75.0f) {
    q_v = 1.025488186204f * q_v + 1.116100545294e-04f;
  }
  return q_v;
}

// ======================= Save to File =======================
void saveToFile(float qv) {
  if (!dataFile) return;
  dataFile.print(values[0], 0); dataFile.print(",");
  dataFile.print(values[1], 0); dataFile.print(",");
  dataFile.print(U0, 3);        dataFile.print(",");
  dataFile.print(values[2], 2); dataFile.print(",");
  dataFile.print(values[3], 2); dataFile.print(",");
  dataFile.println(qv, 6);
  dataFile.flush();
}

// ======================= RUN Command =======================
void sendRun() {
  Serial.print("RUN:");
  Serial.print(values[0], 0); Serial.print(",");
  Serial.print(values[1], 0); Serial.print(",");
  Serial.println(U0, 2); // CR/LF from println
}


// ======================= Handle Key =======================
void handleKey(uint8_t k) {
  // Abort/Reset with key 15
  if (k == 15) {
    sweeping = false;
    currentInput = "";
    for (int i = 0; i < 4; ++i) values[i] = 0;
    stage = ASK_D;
    lcd.clear(); lcd.print("Aborted/Reset");
    delay(800);
    promptUser();
    return;
  }

  switch (k) {
    case 1 ... 9:                      // digits 1..9
      currentInput += char('0' + k);
      break;
    case 10:                           // 0
      currentInput += '0';
      break;
    case 11:                           // decimal point
      if (currentInput.indexOf('.') == -1) currentInput += '.';
      break;
    case 14:                           // backspace
      if (currentInput.length()) currentInput.remove(currentInput.length()-1);
      break;

    case 12: {                         // Enter
      if (!currentInput.length()) break;
      float tempInput = currentInput.toFloat();

      if (stage == ASK_D) {
        if (tempInput != 100 && tempInput != 150) {
          lcd.clear(); lcd.print("Not defined yet");
          delay(1200); promptUser(); return;
        }
        values[ASK_D] = tempInput;
        stage = ASK_d;
        promptUser();
      } else if (stage == ASK_d) {
        // Validate d based on D
        if ((values[ASK_D] == 100 && tempInput != 50) ||
            (values[ASK_D] == 150 && tempInput != 75)) {
          lcd.clear(); lcd.print("Not defined yet");
          delay(1200); promptUser(); return;
        }
        values[ASK_d] = tempInput;

        // Decide sweep params
        if (values[ASK_D] == 100) { U0_start = 0.5f; U0_end = 3.0f; U0_step = 0.1f; }
        else                       { U0_start = 1.0f; U0_end = 3.5f; U0_step = 0.1f; }

        U0 = U0_start;
        sweeping = true;
        sendRun();
      }
      currentInput = "";
      break;
    }
  }

  // Echo current input
  if (k != 12) {
    lcd.setCursor(0,1);
    lcd.print("                ");
    lcd.setCursor(0,1);
    lcd.print(currentInput);
  }
}

// ======================= Setup =======================
void setup() {
  Serial.begin(115200);
  Mon.begin(115200);      // your PC terminal on HW-598 COM

  lcd.init(); lcd.backlight();

  pinMode(KP_SCL_PIN, OUTPUT); digitalWrite(KP_SCL_PIN, HIGH);
  pinMode(KP_SDO_PIN, INPUT);

  if (!SD.begin(SD_CS_PIN)) {
    lcd.clear(); lcd.print("SD init failed!");
    while (1) {}
  }

  askFileMode();        // <<<< choose new/old file here
  promptUser();
}

// ======================= Loop =======================
void loop() {
  // Keypad
  uint8_t k = scanKeypad();
  if (k && k != lastKey) handleKey(k);
  lastKey = k;

  // Read MATLAB reply (P1,P2)
  if (sweeping && Serial.available() > 0) {
    String line = Serial.readStringUntil('\n');  // expects CR/LF from MATLAB writeline
    line.trim();
    int sep = line.indexOf(',');
    if (sep != -1) {
      values[2] = line.substring(0, sep).toFloat();         // P1
      values[3] = line.substring(sep + 1).toFloat();        // P2

      // 1) Show P1 & P2 for 30 seconds
      lcd.clear();
      lcd.setCursor(0,0); lcd.print("P1:"); lcd.print(values[2], 0); lcd.print(" Pa");
      lcd.setCursor(0,1); lcd.print("P2:"); lcd.print(values[3], 0); lcd.print(" Pa");
      delay(30000); // 30 seconds, blocking on purpose per your request

      // 2) Compute q_v and show it (and keep it there)
      float qv = computeQv(values[0], values[1], values[2], values[3]);
      if (!isnan(qv)) {
        saveToFile(qv);
        printBoth("Flow rate: ");
        printBothFloat(qv, 6);
        printlnBoth(" m^3/s");

        lcd.clear();
        lcd.setCursor(0,0); lcd.print("q_v:");
        lcd.setCursor(5,0); lcd.print(qv, 6);
        lcd.setCursor(0,1); lcd.print("m3/s");
      } else {
        printlnBoth("Delta P <= 0. Skipping.");
        lcd.clear();
        lcd.setCursor(0,0); lcd.print("DeltaP <= 0");
        lcd.setCursor(0,1); lcd.print("Check P1,P2");
      }

      // 3) Move to next U0, but DO NOT change the LCD (leave q_v showing)
      U0 += U0_step;
      if (U0 <= U0_end + 1e-6f) {
        delay(50);
        sendRun(); // No LCD writes inside sendRun()
      } else {
        sweeping = false;
        lcd.clear(); lcd.print("Sweep done");
        delay(1200);
        printlnBoth("Sweep done");
        // Reset to next run
        stage = ASK_D;
        currentInput = "";
        for (int i = 0; i < 4; ++i) values[i] = 0;
        promptUser();
      }

    }
  }

  delay(40);
}