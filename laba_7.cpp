#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <random>
#include <string>
#include <sstream>
#include <complex>
#include <algorithm>
#include <functional>

using namespace std;

const double PI = 3.14159265358979323846;

const int n = 9;
const int N = 1 << n; // 512 точек
const double A = 2.94;
const double B = 0.27;
const double w1 = 2.0;     // первая частота
const double w2 = 197.0;   // вторая частота
const double sigma_noise = 197.0; // стандартное отклонение шума

struct Filters {
    vector<double> h;  // низкочастотный фильтр (масштабирующая функция)
    vector<double> g;  // высокочастотный фильтр (вейвлет)
};


void saveTxt(const string& filename, const vector<double>& data) {
    ofstream f(filename);
    if (!f.is_open()) {
        cerr << "Ошибка: не могу открыть файл " << filename << endl;
        return;
    }

    for (size_t i = 0; i < data.size(); i++) {
        f << i << " " << data[i] << endl;
    }
    f.close();
}
// Фильтры Хаара
Filters getHaar() {
    Filters f;
    f.h = { 1.0 / sqrt(2.0), 1.0 / sqrt(2.0) };
    f.g = { 1.0 / sqrt(2.0), -1.0 / sqrt(2.0) };
    return f;
}

// Фильтры Добеши D6
Filters getD6() {
    Filters f;
    f.h = {
        0.3326705529500826,
        0.8068915093110928,
        0.4598775021184915,
       -0.1350110200102546,
       -0.0854412738820267,
        0.0352262918857095
    };

    f.g = {
        0.0352262918857095,
        0.0854412738820267,
       -0.1350110200102546,
       -0.4598775021184915,
        0.8068915093110928,
       -0.3326705529500826
    };
    return f;
}

// Свертка с фильтром и децимация
vector<double> convolveAndDownsample(const vector<double>& signal, const vector<double>& filter) {
    int L = filter.size();
    int N = signal.size();
    int M = N / 2;
    vector<double> result(M, 0.0);

    for (int i = 0; i < M; i++) {
        for (int k = 0; k < L; k++) {
            int idx = 2 * i + k;
            // Периодическое продолжение
            if (idx >= N) idx = idx % N;
            if (idx < 0) idx = N + idx;
            result[i] += filter[k] * signal[idx];
        }
    }

    return result;
}

// Декомпозиция сигнала
void decomposeFIR(const vector<double>& signal, const Filters& f,
    vector<double>& low, vector<double>& high) {
    low = convolveAndDownsample(signal, f.h);
    high = convolveAndDownsample(signal, f.g);
}

// Реконструкция сигнала
vector<double> reconstructFIR(const vector<double>& coeffs, const Filters& f) {
    int M = coeffs.size();
    int N = M * 2;
    vector<double> result(N, 0.0);

    // Интерполяция и свертка
    for (int i = 0; i < M; i++) {
        for (int k = 0; k < f.h.size(); k++) {
            int idx = 2 * i + k;
            if (idx < N) {
                result[idx] += f.h[k] * coeffs[i];
            }
        }
    }

    return result;
}
void processShannon(const vector<double>& signal) {
    int current_size = signal.size();
    vector<double> current = signal;

    // Вычисляем БПФ сигнала
    vector<complex<double>> spectrum(current_size);
    for (int k = 0; k < current_size; k++) {
        complex<double> sum(0, 0);
        for (int n = 0; n < current_size; n++) {
            double angle = -2.0 * PI * k * n / current_size;
            sum += current[n] * complex<double>(cos(angle), sin(angle));
        }
        spectrum[k] = sum;
    }

    // 4-уровневая декомпозиция
    for (int level = 1; level <= 4; level++) {
        int new_size = current_size / 2;
        vector<double> low_freq(new_size, 0.0);
        vector<double> high_freq(new_size, 0.0);

        // Разделяем частоты для текущего уровня
        vector<complex<double>> low_spectrum(current_size, 0.0);
        vector<complex<double>> high_spectrum(current_size, 0.0);

        for (int k = 0; k < current_size; k++) {
            if (k < current_size / 4 || k >= 3 * current_size / 4) {
                low_spectrum[k] = spectrum[k];
            }
            else {
                high_spectrum[k] = spectrum[k];
            }
        }

        // Обратное БПФ и децимация
        for (int i = 0; i < new_size; i++) {
            complex<double> low_sum(0, 0);
            complex<double> high_sum(0, 0);

            for (int k = 0; k < current_size; k++) {
                double angle = 2.0 * PI * k * i / current_size;
                complex<double> factor(cos(angle), sin(angle));
                low_sum += low_spectrum[k] * factor;
                high_sum += high_spectrum[k] * factor;
            }

            low_freq[i] = low_sum.real() / current_size;
            high_freq[i] = high_sum.real() / current_size;
        }

        // Сохраняем результаты
        saveTxt("shannon_low_" + to_string(level) + ".txt", low_freq);
        saveTxt("shannon_high_" + to_string(level) + ".txt", high_freq);

        // Готовимся к следующему уровню
        current = low_freq;
        current_size = new_size;

        // Пересчитываем спектр для следующего уровня
        spectrum.resize(current_size);
        for (int k = 0; k < current_size; k++) {
            complex<double> sum(0, 0);
            for (int n = 0; n < current_size; n++) {
                double angle = -2.0 * PI * k * n / current_size;
                sum += current[n] * complex<double>(cos(angle), sin(angle));
            }
            spectrum[k] = sum;
        }

        // Реконструкция для текущего уровня
        vector<double> reconstructed = current;
        for (int r = 0; r < level; r++) {
            // Интерполяция в 2 раза
            vector<double> interpolated(reconstructed.size() * 2, 0.0);
            for (size_t j = 0; j < reconstructed.size(); j++) {
                interpolated[2 * j] = reconstructed[j];
            }
            reconstructed = interpolated;
        }
        saveTxt("shannon_p_" + to_string(level) + ".txt", reconstructed);
    }
}

int main() {

    system("COLOR F0");
    setlocale(LC_ALL, "RU");

    vector<double> z(N);

    double A = 2.94, B = 0.27;
    double w1 = 2.0, w2 = 197.0, phi = PI / 2.0;

    bool isTask6 = true;

    default_random_engine gen;
    normal_distribution<double> noise(0.0, sigma_noise);

    if (!isTask6) {
        // Задание 2-5: сигнал только на определенных интервалах
        for (int j = 0; j < N; j++) {
            if ((j >= N / 4 && j <= N / 2) || (j > 3 * N / 4))
                z[j] = A + B * cos(2 * PI * w2 * j / N);
            else
                z[j] = 0;
            z[j] += noise(gen);  // Добавляем шум
        }
    }
    else {
        // Задание 6: сигнал на всем интервале с двумя частотами
        for (int j = 0; j < N; j++) {
            z[j] = A * cos(2.0 * PI * w1 * j / N + phi) + B * cos(2.0 * PI * w2 * j / N);
            z[j] += noise(gen);  // Добавляем шум
        }
    }

    saveTxt("signal.txt", z);


    cout << "4-ЭТАПНЫЙ ВЕЙВЛЕТ-АНАЛИЗ СИГНАЛА" << endl;

    auto runMRA = [&](string name, const Filters& f) {
        cout << "\nВейвлет " << name << ":" << endl;
        vector<double> cur = z;
        for (int j = 1; j <= 4; j++) {
            vector<double> l, h;
            decomposeFIR(cur, f, l, h);

            // Сохраняем с другими именами файлов
            saveTxt(name + "_a" + to_string(j) + ".txt", l);  // аппроксимация
            saveTxt(name + "_d" + to_string(j) + ".txt", h);  // детализация

            // Статистика
            cout << "  Этап " << j << ":" << endl;
            cout << "    Аппроксимация: N=" << l.size()
                << ", min=" << *min_element(l.begin(), l.end())
                << ", max=" << *max_element(l.begin(), l.end()) << endl;
            cout << "    Детализация: N=" << h.size()
                << ", min=" << *min_element(h.begin(), h.end())
                << ", max=" << *max_element(h.begin(), h.end()) << endl;

            // Реконструкция для текущего уровня
            vector<double> p = l;
            for (int r = 0; r < j; r++) {
                p = reconstructFIR(p, f);
                // Обрезаем до нужного размера
                if (p.size() > N) {
                    p.resize(N);
                }
            }
            saveTxt(name + "_p_" + to_string(j) + ".txt", p);
            cur = l;
        }
        };

    // Запуск MRA для разных вейвлетов
    cout << "\nВыполнение вейвлет-анализа..." << endl;

    // Хаар
    runMRA("haar", getHaar());

    // Добеши D6
    runMRA("db6", getD6());

    // Шеннон (особая обработка)
    cout << "\nВейвлет Шеннон:" << endl;
    processShannon(z);

    // Переименовываем файлы Шеннона для единообразия
    for (int j = 1; j <= 4; j++) {
        // Переименовываем low -> a, high -> d
        ifstream fin1("shannon_low_" + to_string(j) + ".txt");
        ifstream fin2("shannon_high_" + to_string(j) + ".txt");

        if (fin1.good() && fin2.good()) {
            fin1.close();
            fin2.close();

            // Переименовываем файлы
            string old_low = "shannon_low_" + to_string(j) + ".txt";
            string old_high = "shannon_high_" + to_string(j) + ".txt";
            string new_a = "shannon_a" + to_string(j) + ".txt";
            string new_d = "shannon_d" + to_string(j) + ".txt";

            // Копируем содержимое
            ifstream src_low(old_low, ios::binary);
            ofstream dst_low(new_a, ios::binary);
            dst_low << src_low.rdbuf();
            src_low.close();
            dst_low.close();

            ifstream src_high(old_high, ios::binary);
            ofstream dst_high(new_d, ios::binary);
            dst_high << src_high.rdbuf();
            src_high.close();
            dst_high.close();

            // Удаляем старые файлы
            remove(old_low.c_str());
            remove(old_high.c_str());

            cout << "  Этап " << j << " переименован" << endl;
        }
    }

    // Сохраняем аппроксимацию 4 уровня отдельно
    vector<double> haar_a4, db6_a4, shannon_a4;

    // Читаем последние аппроксимации
    ifstream f1("haar_a4.txt");
    ifstream f2("db6_a4.txt");
    ifstream f3("shannon_a4.txt");

    if (f1.good() && f2.good() && f3.good()) {
        // Файлы уже созданы
        cout << "  Файлы аппроксимации 4 уровня уже существуют" << endl;
    }
    else {
        // Создаем файлы аппроксимации 4 уровня
        vector<double> data;
        for (int i = 0; i < N / 16; i++) data.push_back(0.0);  // заглушка

        saveTxt("haar_a4.txt", data);
        saveTxt("db6_a4.txt", data);
        saveTxt("shannon_a4.txt", data);
    }

    cout << "Созданы файлы:" << endl;
    cout << "  signal.txt - исходный сигнал" << endl;

    for (string wavelet : {"haar", "db6", "shannon"}) {
        cout << "  " << wavelet << "_d1.txt, " << wavelet << "_d2.txt, "
            << wavelet << "_d3.txt, " << wavelet << "_d4.txt - детализация" << endl;
        cout << "  " << wavelet << "_a1.txt, " << wavelet << "_a2.txt, "
            << wavelet << "_a3.txt, " << wavelet << "_a4.txt - аппроксимация" << endl;
        for (int j = 1; j <= 4; j++) {
            cout << "  " << wavelet << "_p_" << j << ".txt - восстановление уровня " << j << endl;
        }
    }

    cout << (isTask6 ? "\nПункт 6 выполнен!" : "\nПункты 2-5 выполнены!") << endl;
    return 0;

}

