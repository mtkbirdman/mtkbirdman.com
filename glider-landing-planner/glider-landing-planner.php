<?php
/**
 * Plugin Name: Glider Landing Planner
 * Description: 教育・机上検討用グライダー着陸経路シミュレーター。[glider_landing_planner]で固定ページへ表示します。
 * Version: 0.1.2
 * Requires PHP: 7.4
 * Text Domain: glider-landing-planner
 * License: GPL-2.0-or-later
 */
if (!defined('ABSPATH')) { exit; }

function glp_register_assets() {
    $base = plugin_dir_url(__FILE__);
    $ver = '0.1.2';
    wp_register_style('glp-leaflet', 'https://unpkg.com/leaflet@1.9.4/dist/leaflet.css', array(), '1.9.4');
    wp_register_style('glp-app', $base . 'assets/css/app.css', array('glp-leaflet'), $ver);
    wp_register_style('glp-print', $base . 'assets/css/print.css', array('glp-app'), $ver, 'print');
    wp_register_script('glp-leaflet', 'https://unpkg.com/leaflet@1.9.4/dist/leaflet.js', array(), '1.9.4', true);
    wp_register_script('glp-app-module', $base . 'assets/js/main.js', array('glp-leaflet'), $ver, true);
}

function glp_enqueue_assets() {
    glp_register_assets();
    wp_enqueue_style('glp-leaflet');
    wp_enqueue_style('glp-app');
    wp_enqueue_style('glp-print');
    wp_enqueue_script('glp-leaflet');
    wp_enqueue_script('glp-app-module');
}

function glp_early_enqueue() {
    if (!is_singular()) { return; }
    $post = get_post();
    if ($post && isset($post->post_content) && has_shortcode($post->post_content, 'glider_landing_planner')) {
        glp_enqueue_assets();
    }
}
add_action('wp_enqueue_scripts', 'glp_early_enqueue', 15);

function glp_render_shortcode($atts = array()) {
    glp_enqueue_assets();
    return '<div class="glp-app-root"><p>Glider Landing Planner を読み込んでいます…</p></div>';
}
add_shortcode('glider_landing_planner', 'glp_render_shortcode');

function glp_module_script_tag($tag, $handle, $src) {
    if ($handle !== 'glp-app-module') { return $tag; }
    // No inline JavaScript. Native ES modules load relative .js files from /assets/js/.
    return str_replace('<script ', '<script type="module" ', $tag);
}
add_filter('script_loader_tag', 'glp_module_script_tag', 10, 3);
