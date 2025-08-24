#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <print>
#include <string>
#include <vector>

#if defined(_WIN32)
#include <Windows.h>
#endif

#include <smath.hpp>
using smath::Vec2;

enum Color : uint8_t {
  CLR_NONE = 0,  // default
  CLR_AXES = 90, // bright black (gray)
  CLR_A = 32,    // green
  CLR_B = 33,    // yellow
  CLR_P = 35,    // magenta
  CLR_DOT = 36   // cyan
};

struct Cell {
  char ch{' '};
  uint8_t color{CLR_NONE};
  int prio{0};
};

struct Canvas {
  int w, h;
  std::vector<Cell> pix;
  Canvas(int W, int H) : w(W), h(H), pix(W * H) {}

  void put(int x, int y, char c, int prio, uint8_t color) {
    if (x < 0 || x >= w || y < 0 || y >= h)
      return;
    Cell &cell = pix[y * w + x];
    if (prio >= cell.prio) {
      cell.ch = c;
      cell.prio = prio;
      cell.color = color;
    }
  }
  void hline(int y, char c, int prio, uint8_t color) {
    for (int x = 0; x < w; ++x)
      put(x, y, c, prio, color);
  }
  void vline(int x, char c, int prio, uint8_t color) {
    for (int y = 0; y < h; ++y)
      put(x, y, c, prio, color);
  }

  void line_dir(int x0, int y0, int x1, int y1, int prio, uint8_t color) {
    int dx = std::abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
    int dy = -std::abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    int px = x0, py = y0;
    put(px, py, '+', prio, color);

    while (true) {
      if (px == x1 && py == y1)
        break;
      int e2 = 2 * err;
      int nx = px, ny = py;
      if (e2 >= dy) {
        err += dy;
        nx += sx;
      }
      if (e2 <= dx) {
        err += dx;
        ny += sy;
      }

      int sdx = nx - px;
      int sdy = ny - py;
      int adx = std::abs(sdx);
      int ady = std::abs(sdy);
      char g;
      if (adx >= 2 * ady)
        g = '-';
      else if (ady >= 2 * adx)
        g = '|';
      else
        g = ((sdx > 0 && sdy > 0) || (sdx < 0 && sdy < 0)) ? '\\' : '/';

      put(nx, ny, g, prio, color);
      px = nx;
      py = ny;
    }
  }

  void flush() const {
    uint8_t cur = 255;
    for (int y = 0; y < h; ++y) {
      for (int x = 0; x < w; ++x) {
        const Cell &c = pix[y * w + x];
        if (c.color != cur) {
          if (c.color == CLR_NONE)
            std::print("\x1b[0m");
          else
            std::print("\x1b[{}m", c.color);
          cur = c.color;
        }
        std::print("{}", c.ch);
      }
      std::println("");
    }
    std::print("\x1b[0m");
  }
};

struct Mapper {
  int W, H;
  double scale;
  int cx, cy;
  Mapper(int w, int h, double s) : W(w), H(h), scale(s), cx(w / 2), cy(h / 2) {}
  std::pair<int, int> map(Vec2 v) const {
    int x = cx + (int)std::llround(v.x() * scale);
    int y = cy - (int)std::llround(v.y() * scale);
    return {x, y};
  }
};

int main() {
#if defined(_WIN32)
  HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
  if (h != INVALID_HANDLE_VALUE) {
    DWORD m;
    if (GetConsoleMode(h, &m))
      SetConsoleMode(h, m | ENABLE_VIRTUAL_TERMINAL_PROCESSING);
  }
#endif

  std::print("Enter vector a (ax ay): ");
  double ax, ay;
  if (!(std::cin >> ax >> ay))
    return 0;
  std::print("Enter vector b (bx by): ");
  double bx, by;
  if (!(std::cin >> bx >> by))
    return 0;

  Vec2 a{(float)ax, (float)ay};
  Vec2 b{(float)bx, (float)by};
  Vec2 p = a.project_onto(b);

  std::println("\nResults:");
  std::println("a = {}", a);
  std::println("b = {}", b);
  std::println("proj_b(a) = p = {}", p);
  std::println("|a|={}  |b|={}  |p|={}", a.magnitude(), b.magnitude(),
               p.magnitude());
  double dot = a.dot(b);
  double ang = (a.magnitude() > 0 && b.magnitude() > 0)
                   ? std::acos(std::clamp(dot / (a.magnitude() * b.magnitude()),
                                          -1.0, 1.0)) *
                         180.0 / M_PI
                   : 0.0;
  std::println("a·b={}  angle(a,b)={} deg\n", dot, ang);

  const int W = 73, H = 33;
  double maxr = std::max({(double)a.magnitude(), (double)b.magnitude(),
                          (double)p.magnitude(), 1.0});
  double usable = 0.45 * std::min(W, H);
  double scale = usable / maxr;

  Canvas cvs(W, H);
  Mapper mp(W, H, scale);

  int pr_axes = 1;
  auto [cx, cy] = mp.map({0, 0});
  cvs.hline(H / 2, '-', pr_axes, CLR_AXES);
  cvs.vline(W / 2, '|', pr_axes, CLR_AXES);
  cvs.put(W / 2, H / 2, 'O', pr_axes + 1, CLR_AXES);

  auto [ox, oy] = mp.map({0, 0});
  auto [ax1, ay1] = mp.map(a);
  auto [bx1, by1] = mp.map(b);
  auto [px1, py1] = mp.map(p);

  int pr_a = 3, pr_b = 3, pr_p = 4;

  cvs.line_dir(ox, oy, ax1, ay1, pr_a, CLR_A);
  cvs.line_dir(ox, oy, bx1, by1, pr_b, CLR_B);
  cvs.line_dir(ox, oy, px1, py1, pr_p, CLR_P);

  if (!a.approx_equal(p) && b.magnitude() > 0) {
    int steps = std::max(std::abs(ax1 - px1), std::abs(ay1 - py1));
    for (int i = 0; i <= steps; ++i) {
      double t = steps ? double(i) / steps : 0.0;
      int x = (int)std::llround(px1 + t * (ax1 - px1));
      int y = (int)std::llround(py1 + t * (ay1 - py1));
      if (i % 2 == 0)
        cvs.put(x, y, '.', 2, CLR_DOT);
    }
  }

  auto place_label = [&](Vec2 v, const char *txt, uint8_t c, int pr) {
    auto tip = mp.map(v * 1.05f);
    for (int i = 0; txt[i]; ++i)
      cvs.put(tip.first + i, tip.second, txt[i], pr, c);
  };
  place_label(a, "A", CLR_A, pr_a + 1);
  place_label(b, "B", CLR_B, pr_b + 1);
  place_label(p, "P", CLR_P, pr_p + 1);

  std::println("Legend: axes(gray), a(green), b(yellow), p=proj(magenta), "
               "dotted=cross-perp\n");
  cvs.flush();
  return 0;
}
