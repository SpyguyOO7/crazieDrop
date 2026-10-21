/**
 * @file main.c
 * @brief Single-shot 10ms inductive coil driver for RP2040.
 *
 * Polls an active-low pushbutton on GPIO 6 (with external 500-ohm pull-up)
 * and drives an NMOS gate on GPIO 5 high for 10ms per press. Waiting for
 * button release before re-arming prevents continuous coil energization
 * if the button is held down.
 */

#include <Arduino.h>
#include <stdbool.h>
#include <stdint.h>
#include "hardware/gpio.h"
#include "pico/time.h"

#define COIL_PIN          5U
#define BUTTON_PIN        6U
#define PULSE_DURATION_MS 100U
#define DEBOUNCE_DELAY_MS 25U
#define REFIRE_DELAY_MS   500U

/* Stringize helper so the pulse duration macro can be printed */
#define __STR(x) #x
#define __EXPAND_STR(x) __STR(x)

/**
 * @brief Configures GPIO directions and ensures the coil starts de-energized.
 */
static void init_hardware(void) {
    gpio_init(COIL_PIN);
    gpio_set_dir(COIL_PIN, GPIO_OUT);
    gpio_put(COIL_PIN, false);

    gpio_init(BUTTON_PIN);
    gpio_set_dir(BUTTON_PIN, GPIO_IN);
    gpio_disable_pulls(BUTTON_PIN); /* External 500-ohm pull-up is used */
}

/**
 * @brief Drives the NMOS gate high for a fixed duration, then pulls it low.
 */
static void pulse_coil(void) {
    Serial.println("Button pressed, sending "
                   __EXPAND_STR(PULSE_DURATION_MS) " ms pulse");
    gpio_put(COIL_PIN, true);
    sleep_ms(PULSE_DURATION_MS);
    gpio_put(COIL_PIN, false);
}

/**
 * @brief Waits for a debounced press-and-release cycle to trigger one pulse.
 *
 * Since this microcontroller performs no other concurrent tasks, blocking
 * delays during the 10ms pulse and 25ms debounce windows eliminate state
 * tracking variables while guaranteeing single-shot behavior.
 */
static void poll_button_once(void) {
    /* Button is active-low; exit early if unpressed (HIGH) */
    if (gpio_get(BUTTON_PIN)) {
        return;
    }

    /* Confirm the press was not transient contact bounce */
    sleep_ms(DEBOUNCE_DELAY_MS);
    if (gpio_get(BUTTON_PIN)) {
        return;
    }

    pulse_coil();

    /* Block until the user releases the button, then debounce the release */
    while (!gpio_get(BUTTON_PIN)) {
        tight_loop_contents();
    }
    sleep_ms(DEBOUNCE_DELAY_MS);

    /* Enforce minimum time between pulses */
    sleep_ms(REFIRE_DELAY_MS);
}

void setup(void) {
    Serial.begin(115200);
    init_hardware();
}

void loop(void) {
    poll_button_once();
}