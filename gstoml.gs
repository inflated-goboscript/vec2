%define LOAD_GSTOML_FP(fp) list gstoml = fp
%define LOAD_GSTOML LOAD_GSTOML_FP(file ```goboscript.toml```)

onflag {
    _gstoml_parse_toml;
}

struct TWConfig {
    frame_rate=30
}

proc _gstoml_parse_toml {

}
