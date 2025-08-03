# vec2.gs

> 2d linear algebra

This is a 2D vector library which is built for [goboscript](https://github.com/aspizu/goboscript).
It is designed to be used with [inflator](https://github.com/faretek1/inflator).

Implements
- A vec2 struct
- A 2x2 matrix struct

## Credits

- 3blue1brown's youtube series on Linear algebra

## Installation

Make sure you have inflator installed

`inflate install https://github.com/FAReTek1/vec2`

add vec2 to your `inflator.toml` config:
```toml
[dependencies]
# ...
vec2 = "https://github.com/FAReTek1/vec2"
```

## Development

use `inflate install -e .`:

1. clone the respository: `git clone https://github.com/FAReTek1/vec2`
2. `cd vec2`
3. `inflate install -e .`
4. `cd test`
5. `inflate`
6. `goboscript build`
7. open `test.sb3`
